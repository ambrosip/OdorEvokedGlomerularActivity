#pragma once
#include <fstream>
#include <mutex>

/// Vidrio represent
namespace vidrio {
/// Platform depenencies
namespace platform {
/// Code using the standard c++ library.
namespace standard {

    using namespace std;

    /// Implementation of the functions requires for io interface.
    /// \see vidrio::platform::io
    struct io_impl {
        io_impl() {} ///< Does nothing. 
                
        /// Opens the file.
        void open(const char* filename) {
            _h.open(filename,fstream::binary);
        }
                
        ///< Closes the file.
        void close() noexcept {
            _h.close();
        }

        /// Reads nbytes into data from the specified offset from the beginning of the file.
        ///
        /// \param[out] data    Buffer for holding the read data.
        /// \param[in]  offset  Location to start reading in the file.
        /// \param[in]  nbytes  The number of bytes to read.
        void read(char* data,uint64_t offset,uint64_t nbytes) {
            while(nbytes){                
                // streamsize might be smaller than uint64_t
                // so iteratively read the potentially smaller 
                // chunks until nbytes are read.
                const auto n=static_cast<streamsize>(nbytes);
                {
                    lock_guard<decltype(_lock)> lock(_lock);
                    _h.seekg(offset).read(data,n);
                }

                offset+=n;
                data+=n;
                nbytes-=n;
            }
        }

        /// \returns the number of bytes in the file.
        uint64_t nbytes() {
            auto pos=_h.tellg();
            _h.seekg(0,ios_base::end);
            auto n=_h.tellg();
            _h.seekg(pos);
            return n;
        }
    private:
        ifstream _h;  ///< file handle
        mutex    _lock; ///< for atomically performing a seek/read pair
    };
} // vidrio::platform::standard

/// Using the pread interface for file io.
namespace linux_pread {

    #include <unistd.h>
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <fcntl.h>

    /// Implementation of the functions requires for io interface.
    /// \see vidrio::platform::io
    struct io_impl {
        io_impl() {} ///< Does nothing. 
                
        /// Opens the file.
        void open(const char* filename) {
#ifdef __APPLE__
            _h=linux_pread::open(filename,O_RDONLY);
#else
            _h=linux_pread::open(filename,O_RDONLY|O_LARGEFILE);
#endif
            if(_h<0)
                throw runtime_error(strerror(errno));
        }
                
        ///< Closes the file.
        void close() noexcept {
            linux_pread::close(_h);
        }

        /// Reads nbytes into data from the specified offset from the beginning of the file.
        ///
        /// \param[out] data    Buffer for holding the read data.
        /// \param[in]  offset  Location to start reading in the file.
        /// \param[in]  nbytes  The number of bytes to read.
        void read(char* data,uint64_t offset,uint64_t nbytes) {
            auto ecode=pread(_h,data,nbytes,offset);
            if(ecode<0)
                throw runtime_error(strerror(errno));

        }

        /// \returns the number of bytes in the file.
        uint64_t nbytes() {
            struct stat stat;
            fstat(_h,&stat);
            return stat.st_size;
        }
    private:
        int _h;  ///< file descriptor 
    };
} //vidrio::platform::linux_vectorized

} // vidrio::platform
} // vidrio
