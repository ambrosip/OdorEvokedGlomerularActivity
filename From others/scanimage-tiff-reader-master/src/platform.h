#pragma once
#include <cstdint>

/// Vidrio represent
namespace vidrio {  
/// Platform depenencies
namespace platform {
        using namespace std;
        
        /// Required io interface.
        ///
        /// Ideally this would describe the minimal depenencies required for the
        /// stuff we want to do in this utility.  It's supposed to be a way of
        /// declaring/injecting depencies.  Not sure that I'm really doing that 
        /// right.
        ///
        /// \tparam T The type that implements required functions.
        ///           For an example see \ref vidrio::platform::win32::io_impl.
        ///           Required methods:
        ///
        ///                 void open(const char* filename);
        ///                 void close() noexcept;
        ///                 void read(void* data,uint64_t offset,unsigned nbytes) const;
        ///                 uint64_t nbytes() const;
        template<typename T>
        struct io {
            explicit io(const char* filename){ t.open(filename); } ///< Opens the file.
            ~io() { t.close(); }                          ///< Closes the file (RAII style).

            /// Reads nbytes into data from the specified offset from the beginning of the file.
            ///
            /// \param[out] data    Buffer for holding the read data.
            /// \param[in]  offset  Location to start reading in the file.
            /// \param[in]  nbytes  The number of bytes to read.
            void read(void* data,uint64_t offset,uint64_t nbytes)  { 
                t.read((char*)data,offset,nbytes); 
            }

            /// \returns the number of bytes in the file.
            uint64_t nbytes() { return t.nbytes(); }
        private:
            T t;
        };

} // vidrio::platform
} // vidrio
