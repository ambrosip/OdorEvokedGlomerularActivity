#pragma once
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

/// Vidrio represent
namespace vidrio {
/// Platform dependencies
namespace platform {
/// Windows specific code
namespace win32 {

    /// Implementation of the functions requires for io interface.
    /// \see vidrio::platform::io
    struct io_impl {
        io_impl() : _h(INVALID_HANDLE_VALUE) {} ///< Does nothing. 
                                                /// Opens the file.
        void open(const char* filename) {
            _h=CreateFile(filename,GENERIC_READ,FILE_SHARE_READ|FILE_SHARE_WRITE,0,
                OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL|FILE_FLAG_OVERLAPPED,0);
            if(_h==INVALID_HANDLE_VALUE)
                throw std::exception("Could not open file.");
        };
                
        ///< Closes the file.
        void close()  {
            CloseHandle(_h);
            _h=INVALID_HANDLE_VALUE;
        }

        /// Reads nbytes into data from the specified offset from the beginning of the file.
        ///
        /// \param[out] data    Buffer for holding the read data.
        /// \param[in]  offset  Location to start reading in the file.
        /// \param[in]  nbytes  The number of bytes to read.
        void read(char* data,uint64_t offset,uint64_t nbytes) const {
            DWORD nread=0;
            OVERLAPPED o={0};
            o.hEvent =CreateEvent(0,TRUE,FALSE,0);
            while(nbytes){
                const auto n=static_cast<DWORD>(nbytes);
                o.Pointer = (PVOID)offset;
                // ReadFile only reads up to 2GB at a time (largest DWORD).
                // So we read in chunks.
                //
                // The use of OVERLAPPED seems a little silly here, since
                // we immediately block, but it does a couple things for us:
                //
                // 1. We opened the file for async use so that this read function
                //    could be called from seperate threads.  This requires use of
                //    OVERLAPPED object.
                // 2. The OVERLAPPED object allows us to atomically read from a 
                //    given offset similar to the linux pread() function.
                ReadFile(_h,data,n,&nread,&o);
                WaitForSingleObject(o.hEvent, INFINITE); // TODO: check return value

                offset+=n;
                data+=n;
                nbytes-=n;
            }
            CloseHandle(o.hEvent);
        }

        /// \returns the number of bytes in the file.
        uint64_t nbytes() const {
            LARGE_INTEGER sz;
            GetFileSizeEx(_h,&sz);
            return sz.QuadPart;
        }
    private:
        HANDLE _h;  ///< file handle        
    };

} // vidrio::platform::win32
} // vidrio::platform
} // vidrio