/// \file Implements the C interface

#include "tiff.reader.h"
#include "platform.h"
#include "tiff.reader.api.h"
using namespace vidrio;
using namespace vidrio::scanimage;
using namespace vidrio::platform;

static const char* AUTHOR="Vidrio Technologies <support@vidriotech.com>";

#undef  assert
#define assert(expr) do{if(!(expr)) throw runtime_error("Assertion failed: " #expr);}while(0)

#ifdef _MSC_VER
#define strdup _strdup
#define snprintf _snprintf
#endif

// Configure io implementation to use.
// TODO: figure out the right way to extern this so I can resolve the platform
//       stuff at link time
// \see vidrio::platform::io
#if defined(_MSC_VER)
    // about a factor of 2 faster to use native file io rather than
    // the standard library (Release build, Win 10 x64, March 2016)
    #include "win32/io.h"
    using reader = tiff::reader<io<win32::io_impl>>;
#else
    #include "linux/io.h"
    typedef tiff::reader<io<linux_pread::io_impl>> reader;
#endif

#if defined(_MSC_VER)
    #define EXPORT __declspec(dllexport)
#endif

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)

/// Identifyies the version of this API.  Use this for support and for
/// adapting your code if the API changes significantly in the future.
/// \return a null-terminated string
EXPORT
const char* ScanImageTiffReader_APIVersion() {
    static int inited=0;
    static char version[1024]={0};
    if(!inited){
        snprintf(version,sizeof(version),
                  "Version %s-%s by %s",STR(GIT_TAG),STR(GIT_HASH),AUTHOR);
        inited=1;
    }
    return version;
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return a ScanImageTiffReader object
EXPORT
ScanImageTiffReader ScanImageTiffReader_Open(const char *filename) {
    try{
        auto p=new reader(filename);
        return ScanImageTiffReader{(void*)p,nullptr};
    } catch(const exception &e) {
        return ScanImageTiffReader{0,strdup(e.what())};
    }
}

EXPORT
void ScanImageTiffReader_Close(ScanImageTiffReader *r) {
    delete (reader*)r->handle;
    r->handle=nullptr;
}


/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 0 on failure, otherwise the number of bytes required to 
///         store the image description string of frame i.
EXPORT
size_t ScanImageTiffReader_GetImageDescriptionSizeBytes(ScanImageTiffReader *r, int i) {    
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->bytesof_header(i);
    }catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 0 on failure, otherwise the number of indexible frames for
///         querying image descriptions.
EXPORT int ScanImageTiffReader_GetImageDescriptionCount(ScanImageTiffReader *r) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->countof_headers();
    }catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return on success, otherwise 0
EXPORT
int ScanImageTiffReader_GetImageDescription(ScanImageTiffReader *r, int i, char *buf, size_t bytesof_buf) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        self->header(i,buf,bytesof_buf);
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
    return 1;
}


/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 0 on failure, otherwise the number of bytes required to 
///         store the concatenated image descriptions from all frames.
EXPORT
size_t ScanImageTiffReader_GetAllImageDescriptionsSizeBytes(ScanImageTiffReader *r) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->bytesof_headers();
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 1 on success, otherwise 0
EXPORT
int ScanImageTiffReader_GetAllImageDescriptions(ScanImageTiffReader *r, char *buf, size_t bytesof_buf) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        self->headers(buf,bytesof_buf);
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
    return 1;
}

/// \returns 0 on failure, otherwise the size of the metadata section in bytes.
///
/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
EXPORT
size_t ScanImageTiffReader_GetMetadataSizeBytes(ScanImageTiffReader *r) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        string m=self->metadata();
        return m.size();
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// \returns 0 on failure, otherwise 1
///
/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string
EXPORT
int ScanImageTiffReader_GetMetadata(ScanImageTiffReader *r, char *buf, size_t bytesof_buf) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        string m=self->metadata();
        assert(m.size()<=bytesof_buf);
        memcpy(buf,m.c_str(),m.size());
        return 1;
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// \returns 0 on failure, otherwise the size of the metadata section in bytes.
///
/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string
EXPORT
size_t ScanImageTiffReader_GetDataSizeBytes(ScanImageTiffReader *r) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->bytesof_data();
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// \returns 0 on failure, otherwise 1
///
/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string
EXPORT
int ScanImageTiffReader_GetData(ScanImageTiffReader *r, char *buf, size_t bytesof_buf) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        self->data(buf,bytesof_buf);
    } catch(const exception &e) {
        r->log=strdup(e.what());
    }
    return r->log==NULL;
}

/// \returns 0 on failure, otherwise a zero'd out struct nd object.
///
/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string
EXPORT
struct nd ScanImageTiffReader_GetShape(ScanImageTiffReader *r) {
    const nd empty={0};
    const reader *self=(reader*)r->handle;
    if(r->log) return empty;
    try{
        return self->shape();
    } catch(const exception &e) {        
        r->log=strdup(e.what());
        return empty;
    }
}

// v1.4 API

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 0 on failure, otherwise the number of indexible frames.
/// \see ScanImageTiffReader_GetImageDescriptionCount
///
/// The "frame count" and the "image description count" are identical.
EXPORT
unsigned ScanImageTiffReader_GetFrameCount(ScanImageTiffReader *r) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->countof_headers();
    }catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
/// \return 0 on failure, otherwise the number of bytes required to store
///         the raw image data in the interval of frames from "beg" up to, 
///         but not including, "end".
EXPORT
size_t ScanImageTiffReader_GetFrameIntervalSizeBytes(ScanImageTiffReader *r, unsigned beg, unsigned end) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        return self->interval_bytesof_data(beg,end);
    } catch(const exception &e) {
        r->log=strdup(e.what());
        return 0;
    }
}

/// On failure, the log field in the returned object will be non-null.
/// It will point to an error string.
///
/// Fills "buf" with raw image data from the interval of frames starting
/// at "beg" and up to, but not including, "end".
///
/// \return 0 on failure, otherwise 1
EXPORT
int ScanImageTiffReader_GetFrameInterval(ScanImageTiffReader *r, unsigned beg, unsigned end, char* buf, size_t bytesof_buf) {
    reader *self=(reader*)r->handle;
    if(r->log) return 0;
    try{
        self->interval(beg,end,buf,bytesof_buf);
    } catch(const exception &e) {
        r->log=strdup(e.what());
    }
    return r->log==NULL;
}
