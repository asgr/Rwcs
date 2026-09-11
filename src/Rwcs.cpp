// While including system headrs, avoid clashes with the old wcsset POXIS
// function that is still exposed in Windows headers.
//
// This old function is *not* exposed if one defines the NO_OLDNAMES
// (MinGW-64bit) or _NO_OLDNAMES (MinGW-32bit) macros. In fact, that's what we
// do for the compilation of wcslib itself. However, these macros cause other
// names to be hidden, and make the C++ system headers internally incompatible,
// turning this into an unfeasible option for compiling this module.
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)
#define wcsset wcsset_
#endif

#include <Rcpp.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

// End the hack above
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)
#undef wcsset
#endif

#include <wcslib.h>
using namespace Rcpp;

// The number of coordinate axes every projection in this file uses. Both
// wcsp2s() and wcss2p() are called with nelem set to this value and with
// buffers sized ncoord * naxis, so the wcsprm handed to them must describe
// exactly naxis axes. See HeaderWcs below for how that is guaranteed.
static const int naxis = 2;

static void enable_wcsperr()
{
  wcserr_enable(1);
  wcsprintf_set(nullptr);
}

// The message wcslib has accumulated in its printf buffer.  enable_wcsperr()
// resets the buffer at the start of every exported call, so read it as-is.  Do
// not reset before reading: wcsprintf_set(nullptr) terminates the buffer at its
// first character, which would discard the message being reported.
static std::string wcs_detail()
{
  const char *buf = wcsprintf_buf();
  return (buf != nullptr && *buf != '\0') ? std::string(buf) : std::string();
}

static void wcs_stop(const std::string &context, int status,
                     const char *const *messages, int nmessages)
{
  std::string msg = "Rwcs: " + context;
  if (status > 0) {
    msg += " (status " + std::to_string(status) + ")";
    if (messages != nullptr && status < nmessages) {
      msg += ": ";
      msg += messages[status];
      msg += ".";
    }
  }
  const std::string detail = wcs_detail();
  if (!detail.empty()) {
    msg += "\n" + detail;
  }
  Rcpp::stop("%s", msg.c_str());
}

static void wcs_stop(const std::string &context)
{
  Rcpp::stop("Rwcs: %s", context.c_str());
}

// Fail if wcsset() has not already recorded a message, since that message is
// far more specific than the generic status string.
static void wcs_stop_set(const std::string &context, const struct wcsprm *wcs,
                         int status)
{
  wcsperr(wcs, "");
  wcs_stop(context, status, wcs_errmsg, 15);
}

static SEXP _wcss2p(struct wcsprm *wcs, NumericVector RA, NumericVector Dec)
{
  const int ncoord = RA.length();
  if (ncoord == 0) {
    return NumericMatrix(0, naxis);
  }

  NumericMatrix world(naxis, ncoord);
  for (int i = 0; i < ncoord; i++) {
    world(0, i) = RA[i];
    world(1, i) = Dec[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix pixel_matrix(naxis, ncoord);
  auto status = wcss2p(wcs, ncoord, naxis,
                       &(world[0]), &(phi[0]), &(theta[0]), &(img[0]),
                       &(pixel_matrix[0]), &(stat[0]));

  if (status) {
    wcsperr(wcs, "");
    // wcss2p() reports a per-coordinate failure as WCSERR_BAD_WORLD, keeps
    // going, and leaves the axis bitmask in stat[] for the caller to recover
    // from.  Any other status means the whole call aborted at wcs.c's cleanup,
    // so stat[] still holds the zeros Rcpp initialises it to, which the caller
    // reads back as "every point projected fine".  There is nothing per-point
    // left to salvage, so raise instead.
    if (status != WCSERR_BAD_WORLD) {
      wcs_stop("world-to-pixel projection failed", status, wcs_errmsg, 15);
    }
    Rcerr << "Failed s2p conversion :(:\n" << wcsprintf_buf();
    return stat;
  }
  return transpose(pixel_matrix);
}

static SEXP _wcsp2s(struct wcsprm *wcs, NumericVector x, NumericVector y)
{
  const int ncoord = x.length();
  if (ncoord == 0) {
    return NumericMatrix(0, naxis);
  }

  NumericMatrix pixel(naxis, ncoord);
  for (int i = 0; i < ncoord; i++) {
    pixel(0, i) = x[i];
    pixel(1, i) = y[i];
  }
  NumericVector phi(ncoord);
  NumericVector theta(ncoord);
  NumericMatrix img(naxis, ncoord);
  IntegerVector stat(ncoord);
  NumericMatrix world_matrix(naxis, ncoord);

  auto status = wcsp2s(wcs, ncoord, naxis,
                       &(pixel[0]), &(img[0]), &(phi[0]), &(theta[0]),
                       &(world_matrix[0]), &(stat[0]));

  if (status) {
    wcsperr(wcs, "");
    // As in _wcss2p(): WCSERR_BAD_PIX is the per-coordinate status and stat[]
    // already holds the axis bitmask, but any other status aborts the call with
    // stat[] still all zeros.
    if (status != WCSERR_BAD_PIX) {
      wcs_stop("pixel-to-world projection failed", status, wcs_errmsg, 15);
    }
    Rcerr << "Failed p2s conversion :(:\n" << wcsprintf_buf();
    return stat;
  }
  return transpose(world_matrix);
}

static void _wcsset(struct wcsprm* wcs,
                    Rcpp::String CTYPE1, Rcpp::String CTYPE2,
                    double CRVAL1, double CRVAL2, double CRPIX1, double CRPIX2,
                    double CD1_1, double CD1_2, double CD2_1, double CD2_2,
                    Rcpp::String RADESYS, int EQUINOX,
                    double PV1_0, double PV1_1, double PV1_2, double PV1_3, double PV1_4,
                    // double PV1_5, double PV1_6, double PV1_7, double PV1_8, double PV1_9, double PV1_10,
                    double PV2_0, double PV2_1, double PV2_2, double PV2_3, double PV2_4, double PV2_5
                    // double PV2_6, double PV2_7, double PV2_8, double PV2_9, double PV2_10
                    )
                    
{
  //setup wcs
  wcs->flag = -1;
  int status = wcsini(1, naxis, wcs);
  if (status) {
    wcs_stop("wcsini() failed while building the WCS from its keyvalues",
             status, wcs_errmsg, 15);
  }

  //insert wcs val
  wcs->crval[0] = CRVAL1;
  wcs->crval[1] = CRVAL2;

  //insert wcs pix
  wcs->crpix[0] = CRPIX1;
  wcs->crpix[1] = CRPIX2;

  //insert wcs cd matrix
  #ifdef HAVE_CD_MATRIX
  wcs->cd[0] = CD1_1;
  wcs->cd[1] = CD1_2;
  wcs->cd[2] = CD2_1;
  wcs->cd[3] = CD2_2;
  #else
  wcs->pc[0] = CD1_1;
  wcs->pc[1] = CD1_2;
  wcs->pc[2] = CD2_1;
  wcs->pc[3] = CD2_2;
  #endif

  //insert ctype safely
  strncpy(wcs->ctype[0], CTYPE1.get_cstring(), sizeof(wcs->ctype[0]) - 1);
  wcs->ctype[0][sizeof(wcs->ctype[0]) - 1] = '\0';
  strncpy(wcs->ctype[1], CTYPE2.get_cstring(), sizeof(wcs->ctype[1]) - 1);
  wcs->ctype[1][sizeof(wcs->ctype[1]) - 1] = '\0';

  //insert radesys and equinox safely
  strncpy(wcs->radesys, RADESYS.get_cstring(), sizeof(wcs->radesys) - 1);
  wcs->radesys[sizeof(wcs->radesys) - 1] = '\0';
  wcs->equinox = EQUINOX;

  //insert wcs pv
#define FILL_PV(WHICH, I, M) \
  if (!R_IsNA(WHICH)) {               \
    wcs->pv[wcs->npv].i = I;          \
    wcs->pv[wcs->npv].m = M;          \
    wcs->pv[wcs->npv].value = WHICH;  \
    wcs->npv++;                       \
  }

  wcs->npv = 0;
  FILL_PV(PV1_0, 1, 0);
  FILL_PV(PV1_1, 1, 1);
  FILL_PV(PV1_2, 1, 2);
  FILL_PV(PV1_3, 1, 3);
  FILL_PV(PV1_4, 1, 4);
  // FILL_PV(PV1_5, 1, 5);
  // FILL_PV(PV1_6, 1, 6);
  // FILL_PV(PV1_7, 1, 7);
  // FILL_PV(PV1_8, 1, 8);
  // FILL_PV(PV1_9, 1, 9);
  // FILL_PV(PV1_10, 1, 10);
  FILL_PV(PV2_0, 2, 0);
  FILL_PV(PV2_1, 2, 1);
  FILL_PV(PV2_2, 2, 2);
  FILL_PV(PV2_3, 2, 3);
  FILL_PV(PV2_4, 2, 4);
  FILL_PV(PV2_5, 2, 5);
  // FILL_PV(PV2_6, 2, 6);
  // FILL_PV(PV2_7, 2, 7);
  // FILL_PV(PV2_8, 2, 8);
  // FILL_PV(PV2_9, 2, 9);
  // FILL_PV(PV2_10, 2, 10);

  //set long and lat axis positions
  wcs->lng = 0;
  wcs->lat = 1;

  status = wcsset(wcs);
  if (status) {
    wcs_stop_set("wcsset() failed while building the WCS from its keyvalues",
                 wcs, status);
  }
}

// The array of wcsprm structs that wcspih() allocates for the several alternate
// coordinate representations in one header. Releasing it in a destructor rather
// than at the call site matters: the constructor of HeaderWcs below can throw
// partway through, and a member object is still destroyed then, while a pair of
// raw data members would leak the array.
class WcsArray {
public:
  WcsArray() = default;
  WcsArray(const WcsArray &) = delete;
  WcsArray &operator=(const WcsArray &) = delete;

  ~WcsArray()
  {
    if (array_ != nullptr) {
      wcsvfree(&nwcs_, &array_);
    }
  }

  int *nwcs() { return &nwcs_; }
  struct wcsprm **array() { return &array_; }

  // The array itself, indexed by the positions wcsidx() reports.
  struct wcsprm *operator[](int i) { return array_ + i; }
  int count() const { return nwcs_; }

private:
  int nwcs_ = 0;
  struct wcsprm *array_ = nullptr;
};

// A two-axis wcsprm carved out of a header with more axes by wcssub(). Only
// freed when it was actually built, since the header may have been two-axis
// already and no copy made.
class WcsSub {
public:
  WcsSub() = default;
  WcsSub(const WcsSub &) = delete;
  WcsSub &operator=(const WcsSub &) = delete;

  ~WcsSub()
  {
    if (active_) {
      wcsfree(&wcs_);
    }
  }

  // Must be called before wcs() is used, and only once.
  void begin()
  {
    wcs_.flag = -1;
    active_ = true;
  }

  struct wcsprm *wcs() { return &wcs_; }

private:
  struct wcsprm wcs_{};
  bool active_ = false;
};

// A WCS read from a FITS header, reduced to the two celestial axes that
// projecting a pixel needs.
//
// wcslib derives the number of coordinate axes from the header itself, taking
// WCSAXESa where present and NAXIS otherwise, so the header of a cube or a 4D
// array yields a wcsprm with three or four axes.  Handing that to wcsp2s() with
// nelem = 2 leaves ncoord and nelem inconsistent with the parsed struct, and
// wcslib's only guard against it is
//
//     ncoord < 1 || (ncoord > 1 && nelem < wcs->naxis)
//
// which a single coordinate slips past because of the ncoord > 1 term.  When it
// does, linp2x() clears naxis doubles into a buffer sized for two: an
// out-of-bounds write that leaves an answer which looks perfectly plausible.
// AddressSanitizer reports it as
//
//     heap-buffer-overflow ... WRITE of size 24 ... in linp2x lin.c:806
//
// for a three-axis header.  Reducing the struct to its celestial pair rather
// than bypassing the check makes the buffer arithmetic correct by construction.
class HeaderWcs {
public:
  HeaderWcs(const Rcpp::String &header, int nkey, int WCSref, int ctrl)
  {
    // wcspih() edits the keyrecords it parses.  Its FLUSH rule compacts the
    // accepted records in place with strncpy(), and a negative ctrl truncates
    // the string at the current position.  Handing it the CHARSXP bytes behind
    // R's shared string pool would corrupt every other copy of that string in
    // the session, so parse a private copy.
    const char *src = header.get_cstring();
    buf_.assign(src, src + strlen(src) + 1);

    int nreject = 0;
    int status = wcspih(buf_.data(), nkey, WCSHDR_all, ctrl,
                        &nreject, array_.nwcs(), array_.array());
    if (status) {
      wcs_stop("failed to read the WCS keyrecords from the supplied header",
               status, wcshdr_errmsg, 6);
    }

    if (array_.count() == 0) {
      wcs_stop("the supplied header contains no WCS keyrecords; supply a "
               "header with CTYPE/CRVAL/CRPIX keys, or use the keyvalues "
               "argument");
    }

    if (WCSref < 0 || WCSref > 26) {
      wcs_stop("the WCS reference must be an index in the range 0-26");
    }

    int alts[27]{};
    status = wcsidx(array_.count(), array_.array(), alts);
    if (status) {
      wcs_stop("wcsidx() failed to index the alternate coordinate "
               "representations in the header", status, wcs_errmsg, 15);
    }

    if (alts[WCSref] < 0) {
      wcs_stop("the header has no WCS with the requested alternate label "
               "(WCSref = " + std::to_string(WCSref) + ")");
    }

    struct wcsprm *parsed = array_[alts[WCSref]];

    if (parsed->naxis == naxis) {
      // Nothing to reduce.  Set it up here rather than letting wcsp2s() do it
      // lazily, so that a header whose axes are individually fine but cannot be
      // combined (a singular CDi_ja, say) reports as a WCS problem on every
      // path through this file, not as an opaque per-point status code.
      int setstat = wcsset(parsed);
      if (setstat) {
        wcs_stop_set("wcsset() failed on the WCS read from the header",
                     parsed, setstat);
      }
      use_ = parsed;
      return;
    }

    // Extract the longitude/latitude pair.  On return, axes[] holds the
    // 1-relative numbers of the source axes that were selected.
    int nsub = naxis;
    int axes[naxis] = {WCSSUB_LONGITUDE, WCSSUB_LATITUDE};
    sub_.begin();
    status = wcssub(1, parsed, &nsub, axes, sub_.wcs());
    if (status) {
      std::string ctx = "wcssub() failed to extract the celestial pair from "
                      + std::to_string(parsed->naxis) + "-axis WCS";
      wcs_stop(ctx, status, wcs_errmsg, 15);
    }

    if (nsub != naxis) {
      wcs_stop("the header's WCS has " + std::to_string(parsed->naxis)
               + " axes but only " + std::to_string(nsub)
               + " celestial axes could be found; exactly 2 are required");
    }

    // wcssub() checks separability against pc, which is still the identity
    // until wcsset() folds CDi_ja into it, so for the common CDi_ja header
    // that check is vacuous and coupling to a dropped axis would be silently
    // discarded.  Verify it here against whichever matrix actually carries the
    // linear transform.
    const double *m = (parsed->cd != nullptr) ? parsed->cd : parsed->pc;
    const int n = parsed->naxis;
    for (int i = 0; i < naxis; i++) {
      const int kept = axes[i] - 1;
      for (int drop = 0; drop < n; drop++) {
        if (drop == kept) continue;
        if (m[kept * n + drop] != 0.0 || m[drop * n + kept] != 0.0) {
          wcs_stop("axis " + std::to_string(kept + 1) + " of this "
                   + std::to_string(n) + "-axis WCS is coupled to axis "
                   + std::to_string(drop + 1) + ", so RA/Dec cannot be "
                   "separated from the remaining axes.  Supply a 2D header, "
                   "or the keyvalues argument, to project this data.");
        }
      }
    }

    // wcssub() deliberately does not run wcsset() on the result.
    status = wcsset(sub_.wcs());
    if (status) {
      std::string ctx = "wcsset() failed on the celestial pair extracted from "
                      + std::to_string(parsed->naxis) + "-axis WCS";
      wcs_stop_set(ctx, sub_.wcs(), status);
    }

    use_ = sub_.wcs();
  }

  HeaderWcs(const HeaderWcs &) = delete;
  HeaderWcs &operator=(const HeaderWcs &) = delete;

  struct wcsprm *get()
  {
    // Central guard on the invariant every projection here relies on.
    if (use_ == nullptr || use_->naxis != naxis) {
      wcs_stop("internal error: the WCS does not describe exactly "
               + std::to_string(naxis) + " axes");
    }
    return use_;
  }

private:
  std::vector<char> buf_;
  WcsArray array_;
  WcsSub sub_;
  struct wcsprm *use_ = nullptr;
};

// [[Rcpp::export]]
SEXP Cwcs_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_0 = NA_REAL, double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL, double PV1_4 = NA_REAL,
              double PV2_0 = NA_REAL, double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL, double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
              )
{
  enable_wcsperr();
  struct wcsprm wcs;
  _wcsset(&wcs, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX,
    PV1_0, PV1_1, PV1_2, PV1_3, PV1_4,
    // PV1_5, PV1_6, PV1_7, PV1_8, PV1_9, PV1_10,
    PV2_0, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5
    // PV2_6, PV2_7, PV2_8, PV2_9, PV2_10
    );
  auto result = _wcss2p(&wcs, RA, Dec);
  wcsfree(&wcs);
  return result;
}

// [[Rcpp::export]]
SEXP Cwcs_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y,
              Rcpp::String CTYPE1 = "RA---TAN", Rcpp::String CTYPE2 = "DEC--TAN",
              double CRVAL1 = 0, double CRVAL2 = 0,
              double CRPIX1 = 0, double CRPIX2 = 0,
              double CD1_1 = 1, double CD1_2 = 0,
              double CD2_1 = 0, double CD2_2 = 1,
              Rcpp::String RADESYS = "ICRS", int EQUINOX = 2000,
              double PV1_0 = NA_REAL, double PV1_1 = NA_REAL, double PV1_2 = NA_REAL, double PV1_3 = NA_REAL, double PV1_4 = NA_REAL,
              double PV2_0 = NA_REAL, double PV2_1 = NA_REAL, double PV2_2 = NA_REAL, double PV2_3 = NA_REAL, double PV2_4 = NA_REAL, double PV2_5 = NA_REAL
              )
{
  enable_wcsperr();
  struct wcsprm wcs;
  _wcsset(&wcs, CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, RADESYS, EQUINOX,
    PV1_0, PV1_1, PV1_2, PV1_3, PV1_4, 
    //PV1_5, PV1_6, PV1_7, PV1_8, PV1_9, PV1_10,
    PV2_0, PV2_1, PV2_2, PV2_3, PV2_4, PV2_5
    //PV2_6, PV2_7, PV2_8, PV2_9, PV2_10
    );
  auto result = _wcsp2s(&wcs, x, y);
  wcsfree(&wcs);
  return result;
}

// [[Rcpp::export]]
SEXP Cwcs_head_p2s(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::String header, 
                   int nkey, int WCSref=0, int ctrl=2)
{
  enable_wcsperr();
  HeaderWcs wcs(header, nkey, WCSref, ctrl);
  return _wcsp2s(wcs.get(), x, y);
}

// [[Rcpp::export]]
SEXP Cwcs_head_s2p(Rcpp::NumericVector RA, Rcpp::NumericVector Dec, Rcpp::String header, 
                   int nkey, int WCSref=0, int ctrl=2)
{
  enable_wcsperr();
  HeaderWcs wcs(header, nkey, WCSref, ctrl);
  return _wcss2p(wcs.get(), RA, Dec);
}
