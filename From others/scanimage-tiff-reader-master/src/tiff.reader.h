// Copyright 2016 Vidrio Technologies, LLC
#include <map>
#include <list>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <future>
#include <cstddef> // for offsetof
#include <string>
#include <cstring>
#include "nd.h"
#include <sstream>

namespace vidrio {
/// ScanImage utilites
namespace scanimage {
/// ScanImage Tiff utilites
namespace tiff{

    using namespace std;

    struct rational  { uint32_t num, den; };
    struct srational {  int32_t num, den; };

    enum Endian {
        LittleEndian,
        BigEndian
    };

	enum MetadataFormat {
		MetadataFormat_Old,
		MetadataFormat_2016_1,
		MetadataFormat_None
	};

#pragma pack(push,1)
    /// Tiff IFD entry structure
    struct tiff_entry{
        uint16_t tag,    ///< Tag identifier
                 type;   ///< Element type identifier
        uint32_t count,  ///< Element count
                 data;   ///< Data section or pointer.  If data doesn't fit in this section, this is an offset to the data's location in the file.
    };

    /// BigTiff IFD entry structure
    struct bigtiff_entry{
        uint16_t tag,    ///< Tag identifier
                 type;   ///< Element type identifier
        uint64_t count,  ///< Element count
                 data;   ///< Data section or pointer.  If data doesn't fit in this section, this is an offset to the data's location in the file.
    };

	/// metadata header block for SI 2016+ files
	struct metadata_header_2016_1 {
		uint32_t
			magic_number,
			version_id,
			header_section_nbytes,
			roigroup_section_nbytes;
	};
#pragma pack(pop)

    template<Endian SourceEndian,typename T> inline void byteswap(T* v);
    template<> inline void byteswap<LittleEndian>(short* v) {}
    template<> inline void byteswap<LittleEndian>(uint16_t* v) {}
    template<> inline void byteswap<LittleEndian>(int* v) {}
    template<> inline void byteswap<LittleEndian>(uint32_t* v) {}
    template<> inline void byteswap<LittleEndian>(uint64_t* v) {}
    template<> inline void byteswap<LittleEndian>(tiff_entry *e){}
    template<> inline void byteswap<LittleEndian>(bigtiff_entry *e){}
#ifdef _MSC_VER
    // Visual Studio builtins
    template<> inline void byteswap<BigEndian>(short* v)    { *v=_byteswap_ushort(*v);}
    template<> inline void byteswap<BigEndian>(uint16_t* v) { *v=_byteswap_ushort(*v);}
    template<> inline void byteswap<BigEndian>(int* v)      { *v=_byteswap_ulong(*v);}
    template<> inline void byteswap<BigEndian>(uint32_t* v) { *v=_byteswap_ulong(*v);}
    template<> inline void byteswap<BigEndian>(uint64_t* v) { *v=_byteswap_uint64(*v);}
#else
    // GCC builtins
    template<> inline void byteswap<BigEndian>(short* v)    { *v=__builtin_bswap16(*v);}
    template<> inline void byteswap<BigEndian>(uint16_t* v) { *v=__builtin_bswap16(*v);}
    template<> inline void byteswap<BigEndian>(int* v)      { *v=__builtin_bswap32(*v);}
    template<> inline void byteswap<BigEndian>(uint32_t* v) { *v=__builtin_bswap32(*v);}
    template<> inline void byteswap<BigEndian>(uint64_t* v) { *v=__builtin_bswap64(*v);}
#endif
    template<> inline void byteswap<BigEndian>(tiff_entry *e) {
        byteswap<BigEndian>(&e->tag);
        byteswap<BigEndian>(&e->type);
        byteswap<BigEndian>(&e->count);
        byteswap<BigEndian>(&e->data);
    }
    template<> inline void byteswap<BigEndian>(bigtiff_entry *e) {
        byteswap<BigEndian>(&e->tag);
        byteswap<BigEndian>(&e->type);
        byteswap<BigEndian>(&e->count);
        byteswap<BigEndian>(&e->data);
    }



    template<Endian SourceEndian,typename T> inline void byteswapv(vector<T> &v) {}
    // if these are ever slow, thse are good candidates for SIMD
    template<> inline void byteswapv<BigEndian>(vector<uint16_t> &v) {
        std::for_each(v.begin(),v.end(),[](uint16_t &r){
            byteswap<BigEndian>(&r);
        });
    }
    template<> inline void byteswapv<BigEndian>(vector<uint32_t> &v) {
        std::for_each(v.begin(),v.end(),[](uint32_t &r){
            byteswap<BigEndian>(&r);
        });
    }
    template<> inline void byteswapv<BigEndian>(vector<uint64_t> &v) {
        std::for_each(v.begin(),v.end(),[](uint64_t &r){
            byteswap<BigEndian>(&r);
        });
    }

    template<Endian SourceEndian,typename T> inline void byteswapv(char* buf,size_t nbytes) {
        const char* end=buf+nbytes;
        for(char* cur=buf;cur<end;cur+=sizeof(T))
            byteswap<SourceEndian,T>((T*)cur);
    }
    template<> inline void byteswapv<LittleEndian,char>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,int16_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,int32_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,int64_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,uint16_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,uint32_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,uint64_t>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,float>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<LittleEndian,double>(char* buf,size_t nbytes) {}
    template<> inline void byteswapv<BigEndian,char>(char* buf,size_t nbytes) {}

    enum tiff_kind_t{
        TiffKindNormal,
        TiffKindBigTiff
    };

    /// Keep interesting tiff tags here.
    enum TiffTag{
        TIFFTAG_IMAGEWIDTH                =256,   ///< The number of columns in the image, i.e., the number of pixels per scanline.
        TIFFTAG_IMAGELENGTH               =257,   ///< The number of rows (sometimes described as scanlines) in the image
        TIFFTAG_BITSPERSAMPLE             =258,
        TIFFTAG_COMPRESSION               =259,
        TIFFTAG_PHOTOMETRICINTERPRETATION =262,
        TIFFTAG_IMAGEDESCRIPTION          =270,   ///< Image description tag from the Tiff standard.
        TIFFTAG_STRIPOFFSETS              =273,
        TIFFTAG_ORIENTATION               =274,   ///< The orientation of the image with respect to the rows and columns
        TIFFTAG_SAMPLESPERPIXEL           =277,   ///< The number of components per pixel. This number is 3 for RGB images, unless extra samples are present. See the ExtraSamples field for further information
        TIFFTAG_ROWSPERSTRIP              =278,
        TIFFTAG_STRIPBYTECOUNTS           =279,
        TIFFTAG_XRESOLUTION               =282,
        TIFFTAG_YRESOLUTION               =283,
        TIFFTAG_PLANARCONFIG              =284,
        TIFFTAG_RESOLUTIONUNIT            =296,
        TIFFTAG_TILEWIDTH                 =322,
        TIFFTAG_TILEHEIGHT                =323,
        TIFFTAG_TILEOFFSETS               =324,
        TIFFTAG_TILEBYTECOUNTS            =325,
        TIFFTAG_SAMPLEFORMAT              =339,   ///< This field specifies how to interpret each data sample in a pixel
    };

    enum TiffResolutionunit {
        TIFFRESOLUTIONUNIT_RELATIVE   =1,
        TIFFRESOLUTIONUNIT_INCH       =2,
        TIFFRESOLUTIONUNIT_CENTIMETER =3,
    };

    enum TiffSampleFormat {
        TIFFSAMPLEFORMAT_UNSIGNED  =1,
        TIFFSAMPLEFORMAT_SIGNED    =2,
        TIFFSAMPLEFORMAT_FLOATING  =3,
        TIFFSAMPLEFORMAT_UNDEFINED =4   ///< A field value of �undefined� is a statement by the writer that it did not know how to interpret the data samples.  Should probably default to unsigbed.
    };

    enum TiffOrientation {
        cr_topleft     =1, ///< (default) The 0th row represents the visual top of the image, and the 0th column represents the visual left-hand side.
        cr_topright    =2,
        cr_bottomright =3,
        cr_bottomleft  =4,
        rc_topleft     =5, ///< The 0th row represents the visual left-hand side of the image, and the 0th column represents the visual top
        rc_topright    =6,
        rc_bottomright =7,
        rc_bottomleft  =8,
    };

    enum TiffType {
        TIFFTYPE_BYTE=1,
        TIFFTYPE_ASCII,
        TIFFTYPE_SHORT,
        TIFFTYPE_LONG,
        TIFFTYPE_RATIONAL,
        TIFFTYPE_SBYTE,
        TIFFTYPE_UNDEFINED,
        TIFFTYPE_SSHORT,
        TIFFTYPE_SLONG,
        TIFFTYPE_SRATIONAL,
        TIFFTYPE_FLOAT,
        TIFFTYPE_DOUBLE,
		TIFFTYPE_IFD4, // aka LONG, seems to have been in the draft spec used by aware systems
        TIFFTYPE_LONG8=16,
        TIFFTYPE_SLONG8,
        TIFFTYPE_IFD8,
    };

    inline uint64_t sizeof_type(TiffType t) {
        switch (t) {
#define XXX(T,n) case TIFFTYPE_##T: return n;
            XXX(BYTE, 1);
            XXX(ASCII, 1);
            XXX(SHORT, 2);
            XXX(LONG, 4);
            XXX(RATIONAL, 8);
            XXX(SBYTE, 1);
            XXX(UNDEFINED, 0);
            XXX(SSHORT, 2);
            XXX(SLONG, 4);
            XXX(SRATIONAL, 8);
            XXX(FLOAT, 4);
            XXX(DOUBLE, 8);
			XXX(IFD4, 4);
            XXX(LONG8, 8);
            XXX(SLONG8, 8);
            XXX(IFD8, 8);
#undef XXX
        default:
            throw runtime_error("Got unrecognized element type in a tiff tag.");
        }
        return 0;
    }
    inline uint64_t sizeof_type(uint16_t t){ return sizeof_type((TiffType)t); }

    struct strip_t {
        uint64_t offset, nbytes;
    };

    struct data_t {
        enum TiffType type;
        uint64_t nelem,offset;
        inline uint64_t nbytes() const { return nelem*sizeof_type(type); }
    };

    struct index_t{
        uint64_t next;               ///< offset of next ifd in file
        map<TiffTag,data_t> tags;
        vector<strip_t>     data;    ///< offsets to strips/tile image data
        unsigned bits_per_sample;
        inline bool has(TiffTag t)const { return tags.find(t) != tags.end(); }
        /// \returns the offset to the last byte of the last strip
        inline size_t end() const {
            uint64_t mx=0;
            for(const auto& s:data)
                mx=std::max(mx,s.nbytes+s.offset);
            return mx;
        }
    };

    #define countof(e) ((sizeof(e)/sizeof(*(e))))
    static void squeeze(nd& s) {
        auto o = 0;
        for (auto i = 0U; i < s.ndim;++i) {
            if(s.dims[i]==1) {
                o -= 1;
            } else {
                s.dims[i+o] = s.dims[i];
            }
        }
        s.ndim += o;
        // Fill remaining higher dims with ones (like it's a high dimensional plane)
        // This makes it simpler to consume the "shape" structure in some cases.
        for(auto i=s.ndim;i<countof(s.dims);++i)
            s.dims[i]=1;
    }
    #undef countof

    static void restride(nd& s) {
        s.strides[0] = 1;
        for (auto i = 0U; i < s.ndim; ++i)
            s.strides[i + 1] = s.dims[i] * s.strides[i];
    }

	//
	// Metadata processing
	//

	template<MetadataFormat> void read_metadata_nbytes(uint32_t* nbytes, function<void(void*,int64_t,uint64_t)> reader);
	template<MetadataFormat> void read_metadata_rawbytes(vector<char> &buf,uint32_t nbytes, function<void(void*,int64_t,uint64_t)> reader);
	template<MetadataFormat> string process_metadata(char *buf,uint64_t n, function<void(void*,int64_t,uint64_t)> reader){ return string(buf,n); } ///< Deserialize the input data stream using libmx

	template<> void read_metadata_nbytes<MetadataFormat_Old>(uint32_t* nbytes, function<void(void*,int64_t,uint64_t)> reader){
		// reader() offset must be seeking from the end...this is super hokey, sorry
		reader(nbytes, sizeof(uint32_t), sizeof(uint32_t));
	}

	template<> void read_metadata_nbytes<MetadataFormat_2016_1>(uint32_t* nbytes, function<void(void*,int64_t,uint64_t)> reader){
		metadata_header_2016_1 metadata_header={0};
		reader(&metadata_header,16,sizeof(metadata_header));
		*nbytes=metadata_header.header_section_nbytes+metadata_header.roigroup_section_nbytes;
	}

	template<> void read_metadata_rawbytes<MetadataFormat_Old>(vector<char> &buf,uint32_t nbytes, function<void(void*,int64_t,uint64_t)> reader){
		// reader() offset must be seeking from the end...this is super hokey, sorry
		reader(buf.data(),sizeof(nbytes) - nbytes,nbytes);
	}

	template<> void read_metadata_rawbytes<MetadataFormat_2016_1>(vector<char> &buf,uint32_t nbytes, function<void(void*,int64_t,uint64_t)> reader){
		reader(buf.data(),32,nbytes);
	}

	static string remove_enclosing_braces(const string& s) {
		auto
			b=s.find("{"),
			e=s.rfind("}");
		return s.substr(b+1,e-1);
	}

	template<> string process_metadata<MetadataFormat_2016_1>(char *buf,uint64_t n, function<void(void*,int64_t,uint64_t)> reader){
		metadata_header_2016_1 metadata_header={0};
		reader(&metadata_header,16,sizeof(metadata_header));
		if(metadata_header.roigroup_section_nbytes==0)
			return string(buf,n);
		int is_json=buf[0]=='{'; // FIXME: fragile.  would be nice to use the version identifier to distinguish encoding types
		if(!is_json) {
			buf[metadata_header.header_section_nbytes-1]='\n'; // replace null character with newline
			return string(buf,n);
		} else {
			string a(buf), // null termination point is not necessarily the same as the byte length
				b(buf+metadata_header.header_section_nbytes);
			a=remove_enclosing_braces(a);
			b=remove_enclosing_braces(b);
			return string("{")+a+string(",")+b+string("}");
		}
	}

    // define a range type to help with iterating over intervals of frames 

    template <class Iter>
    class range {
        Iter b,e;
    public:
        range(Iter b, Iter e) : b(b), e(e) {}
        Iter begin() { return b; }
        Iter end() { return e; }
    };

    template <class Container>
    range<typename Container::iterator> 
    make_range(Container& c, size_t b, size_t e) {
        return range<typename Container::iterator> (c.begin()+b, c.begin()+e);
    }

    /// Opens a tiff file for the specific goal of extracting
    /// image descriptions and the metadata.
    /// \tparam file_t A concrete realization of the \ref platform::io generic type.
    ///                This realizes a platform-specific implementation of the io
    ///                interface required by this class.
    template<typename file_t>
    class reader {
        file_t _file;                    ///< file handle \see platform::io
        vector<index_t> _index;          ///< Index listing where the ifds (and other data) are located in the file.
        tiff_kind_t _kind;               ///< Identifies the tiff as a BigTiff or normal Tiff
        uint64_t _first_ifd;             ///< Location of the start of the first IFD.
        Endian _endian;                  ///< Endian-ness of input file as indicated by the Tiff Magic Number

        unsigned w, h, nplanes, nchan;
        enum nd_type type;

        size_t end_of_image_block;
    public:

        /// Opens a tiff file for the specific goal of extracting
        /// image descriptions and the metadata.
        ///
        /// \param[in] filename  Path to the Tiff file to read.
        ///                      Both BigTiff and Tiff should be accepted.
        reader(const char* filename):
            _file(filename),
            _kind(TiffKindNormal),
            _first_ifd(0)
        {
            identify();
            index();
        }

        nd shape() const {
            nd s={0};
            s.type = type;
            s.ndim = 4;
            s.dims[0] = nchan;
            s.dims[1] = w;
            s.dims[2] = h;
            s.dims[3] = nplanes;
            squeeze(s);
            restride(s);
            return s;
        }

        uint64_t bytesof_data() {
            uint64_t n = 0;
            for(auto ifd:_index) {
                for (auto strip : ifd.data)
                    n += strip.nbytes;
            }
            return n;
        }

        vector<char> data() {
            vector<char> out(bytesof_data());
            data(out.data(),out.size());
            return out;
        }

        void data(char *buf,size_t nbytes){
            switch(_endian){
                case BigEndian:     data<BigEndian>(buf,nbytes); break;
                case LittleEndian:  data<LittleEndian>(buf,nbytes); break;
                default: throw runtime_error("Wrong.");
            }
        }

        template<Endian SourceEndian>
        void data(char *buf,size_t bytesof_buf) {
            auto nbytes = bytesof_data();
            if (!nbytes)
                throw runtime_error("No image data found.  Seems strange.");
            if( bytesof_buf<nbytes )
                throw runtime_error("Input buffer size wasn't large enough.");

            {
                auto* cur=buf;
                list<future<void>> jobs;
                for (auto& ifd : _index) {
                    for (auto& strip : ifd.data) {
                        jobs.push_back(async(launch::async,
                            &file_t::read, &_file, cur, strip.offset, strip.nbytes));
                        cur += strip.nbytes;
                    }
                }
                for(auto &j:jobs)
                    j.get();
            }

            {// Handle endian
                auto* cur=buf;
                for (auto& ifd : _index) {
                    size_t sz=0;
                    for(auto& strip:ifd.data)
                        sz+=strip.nbytes;
                    switch(ifd.bits_per_sample) {
                    case 8:
                        break; // do nothing
                    case 16: byteswapv<SourceEndian, uint16_t>(cur, sz); break;
                    case 24: break;
                    case 32: byteswapv<SourceEndian, uint32_t>(cur, sz); break;
                    case 64: byteswapv<SourceEndian, uint64_t>(cur, sz); break;
                    default:
                        throw runtime_error("Unexpected number of bits per sample.  Not sure what to do.");
                    }

                    cur+=sz;
                }
            }
        }

        uint64_t interval_bytesof_data(size_t b, size_t e) {
            uint64_t n = 0;
            for(auto ifd:make_range(_index,b,e)) {
                for (auto strip : ifd.data)
                    n += strip.nbytes;
            }
            return n;
        }

        vector<char> interval(unsigned b, unsigned e) {
            vector<char> out(interval_bytesof_data(b,e));
            interval(b,e,out.data(),out.size());
            return out;
        }

        void interval(unsigned b, unsigned e, char *buf,size_t nbytes){
            switch(_endian){
                case BigEndian:     interval<BigEndian>(b,e,buf,nbytes); break;
                case LittleEndian:  interval<LittleEndian>(b,e,buf,nbytes); break;
                default: throw runtime_error("Wrong.");
            }
        }

        template<Endian SourceEndian>
        void interval(unsigned b, unsigned e, char *buf,size_t bytesof_buf) {
            auto nbytes = interval_bytesof_data(b,e);
            if (!nbytes)
                throw runtime_error("No image data found.  Seems strange.");
            if( bytesof_buf<nbytes )
                throw runtime_error("Input buffer size wasn't large enough.");

            {
                auto* cur=buf;
                list<future<void>> jobs;
                for (auto& ifd : make_range(_index,b,e)) {
                    for (auto& strip : ifd.data) {
                        jobs.push_back(async(launch::async,
                            &file_t::read, &_file, cur, strip.offset, strip.nbytes));
                        cur += strip.nbytes;
                    }
                }
                for(auto &j:jobs)
                    j.get();
            }

            {// Handle endian
                auto* cur=buf;
                for (auto& ifd : make_range(_index,b,e)) {
                    size_t sz=0;
                    for(auto& strip:ifd.data)
                        sz+=strip.nbytes;
                    switch(ifd.bits_per_sample) {
                    case 8:
                        break; // do nothing
                    case 16: byteswapv<SourceEndian, uint16_t>(cur, sz); break;
                    case 24: break;
                    case 32: byteswapv<SourceEndian, uint32_t>(cur, sz); break;
                    case 64: byteswapv<SourceEndian, uint64_t>(cur, sz); break;
                    default:
                        throw runtime_error("Unexpected number of bits per sample.  Not sure what to do.");
                    }

                    cur+=sz;
                }
            }
        }

        inline size_t bytesof_header(size_t i) {
            auto& ifd = _index.at(i);
            auto it = ifd.tags.find(TIFFTAG_IMAGEDESCRIPTION);
            return (it == ifd.tags.end()) ? 0 : it->second.nbytes();
        }

        /// \returns the number of entries  in the ifd index.
        inline size_t countof_headers() const {
            return _index.size();
        }

        /// Fills buf with the image description from frame i
        /// If you're reading a bunch of image descriptions it's probably
        /// better to read them all (\see headers()), which tries to
        /// read them in parallel.
        void header(size_t i,char *buf,size_t nbytes) {
            auto& ifd = _index.at(i);
            auto it = ifd.tags.find(TIFFTAG_IMAGEDESCRIPTION);
            if (it != ifd.tags.end()) {
                auto & tag = it->second;
                nbytes=std::min(nbytes,(size_t)tag.nbytes());
                _file.read(buf, tag.offset, nbytes);
            }
        }

        inline size_t bytesof_headers() {
            return accumulate(_index.begin(),_index.end(),0ULL,
                [](uint64_t a,const index_t &b){
                    auto it = b.tags.find(TIFFTAG_IMAGEDESCRIPTION);
                    return a + ((it == b.tags.end()) ? 0 : it->second.nbytes());
                });
        }

        /// \returns a list of image description strings
        vector<string> headers() {
            vector<string> hs;
            const auto nbytes=bytesof_headers();
            if (!nbytes)
                return vector<string>();
            vector<char> buf(nbytes,0);
            headers(buf.data(),nbytes);
            // package things up as strings.
            {
                auto* cur=buf.data();
                for(auto ifd: _index)  {
                    auto it = ifd.tags.find(TIFFTAG_IMAGEDESCRIPTION);
                    if (it != ifd.tags.end()) {
                        auto & tag = it->second;
                        //Don't want to use the begin/end iterators here because there might
                        //be trailing nulls.  Want to copy up to the terminating nulls.
                        string s(cur, strlen(buf.data()));
                        cur += tag.nbytes();
                        s.erase(s.find_last_not_of(' ') + 1); // trim
                        hs.push_back(s);
                    }
                }
            }
            return hs;
        }

        /// Reads all the headers into a buffer.  They're concatenated together.
        void headers(char* buf,size_t nbytes) {
            switch(_endian) {
                case BigEndian:    return headers<BigEndian>(buf,nbytes);
                case LittleEndian: return headers<LittleEndian>(buf, nbytes);
                default: throw runtime_error("Wrong.");
            }
        }

        /// headers()  Generalized over endianness
        template<Endian SourceEndian>
        void headers(char* buf,size_t nbytes) {
            auto* cur=buf;
            list<future<void>> jobs;
            for (auto ifd : _index) {
                auto it = ifd.tags.find(TIFFTAG_IMAGEDESCRIPTION);
                if (it != ifd.tags.end()) {
                    auto & tag = it->second;
                    jobs.push_back(async(launch::async,
                        &file_t::read, &_file, cur, tag.offset, tag.nbytes()));
                    cur += tag.nbytes();
                }
            }
            for(auto &j:jobs)
                j.get();
            // This compiles to a no-op because 1-byte elements don't need a
            // byteswap.  But if we need to support wider-character streams in
            // the future we'd want it to kick in without surprises.
            byteswapv<SourceEndian,char>(buf,nbytes);
        }

        /// \returns the metadata as a string.
        string metadata(){
            switch(_endian){
                case BigEndian:    return metadata<BigEndian>();
                case LittleEndian: return metadata<LittleEndian>();
                default: throw runtime_error("Wrong.");
            }
        }

        /// metadata() Generalized over endianness
        template<Endian SourceEndian>
		string metadata(){
			uint32_t nbytes=0;

			MetadataFormat fmt=detect_metadata_format();
			auto reader=[&](void* p,int64_t o,uint64_t n) { _file.read(p,o,n); };
			auto fromend=[&](void* p,int64_t o,uint64_t n) { _file.read(p,_file.nbytes()-o,n); };

			switch(fmt){
				case MetadataFormat_2016_1: read_metadata_nbytes<MetadataFormat_2016_1>(&nbytes,reader); break;
				case MetadataFormat_Old:    read_metadata_nbytes<MetadataFormat_Old>(&nbytes,fromend); break;
				default: return "";
			}
			byteswap<SourceEndian>(&nbytes);
			if(!nbytes)
				return "";
			vector<char> buf(nbytes, 0);
			switch(fmt){
				case MetadataFormat_2016_1: read_metadata_rawbytes<MetadataFormat_2016_1>(buf,nbytes,reader); break;
				case MetadataFormat_Old:    read_metadata_rawbytes<MetadataFormat_Old>(buf,nbytes,fromend); break;
				default: return "";
			}
			byteswapv<SourceEndian>(buf);
			switch(fmt){
				case MetadataFormat_2016_1: return process_metadata<MetadataFormat_2016_1>(buf.data(),buf.size(),reader);
				case MetadataFormat_Old:    return process_metadata<MetadataFormat_Old>(buf.data(),buf.size(),reader);
				default: return "";
			}
		}


    private:

		MetadataFormat detect_metadata_format() {
			metadata_header_2016_1 metadata_header={0};
			_file.read(&metadata_header,16,sizeof(metadata_header));
			// NOTE: Do I need to worry about endian here?
			//       At the moment, if the endian is wrong then the header just won't be decoded.
			//       Not sure how such a file would be generated.
			if(metadata_header.magic_number==0x07030301ULL) {
				switch(metadata_header.version_id) {
					case 3:  // this is the proper id for MetadataFormat_2016_1
					default: // not sure what other values might have been in old files
						return MetadataFormat_2016_1;
				}
			}

			if(end_of_image_block!=_file.nbytes()) // some data was appended
				return MetadataFormat_Old;
			return MetadataFormat_None;
		}



        /// Looks at the first few bytes of the file to check the file type, endianness,
        /// and to identify the kind of Tiff
        void identify() {
            struct {
                short tiff_magic_number;
                short tiff_subtype;
            } data={0};
            static_assert(sizeof(data)==4,"Tiff file id section must be 4 bytes.");
            _file.read(&data,0,sizeof(data));
            switch(data.tiff_magic_number) {
                case 0x4949: _endian=LittleEndian; break;
                case 0x4d4d: _endian=BigEndian; break;
                default:     throw invalid_argument("Invalid tiff.  Endian-ness check returned an invalid id.");
            }
            switch(data.tiff_subtype) {
                case 0x2A00:
                case 0x002A: _kind=TiffKindNormal; break;
                case 0x2B00:
                case 0x002B: _kind=TiffKindBigTiff; break;
                default:
                    throw invalid_argument("Unrecognized Tiff sub-type");
            }
            switch(_endian) {
                case LittleEndian: read_first_ifd_offset<LittleEndian>(); break;
                case BigEndian:    read_first_ifd_offset<BigEndian>();    break;
            }
        }

        template<Endian SourceEndian>
        void read_first_ifd_offset() {
            switch(_kind) {
                case TiffKindNormal:  _file.read(&_first_ifd,4,4); byteswap<SourceEndian>((uint32_t*)&_first_ifd); break;
                case TiffKindBigTiff: _file.read(&_first_ifd,8,8); byteswap<SourceEndian,uint64_t>(&_first_ifd); break;
            }

        }

        /// Iterates over the ifds and build's an index that tells us where the
        /// interesting parts of the file are.
        void index() {
            switch(_endian){
                case LittleEndian:
                    switch(_kind){
                        case TiffKindNormal:  gen_index<uint16_t,uint32_t,tiff_entry,LittleEndian>();    break;
                        case TiffKindBigTiff: gen_index<uint64_t,uint64_t,bigtiff_entry,LittleEndian>(); break;
                    }
                    break;
                case BigEndian:
                    switch(_kind){
                        case TiffKindNormal:  gen_index<uint16_t,uint32_t,tiff_entry,BigEndian>();    break;
                        case TiffKindBigTiff: gen_index<uint64_t,uint64_t,bigtiff_entry,BigEndian>(); break;
                    }
                    break;
            }
            end_of_image_block=0;
            for(const auto& i:_index)
                end_of_image_block=std::max(end_of_image_block,i.end());
        }

        template<typename entry_t,Endian SourceEndian>
        void readtag(index_t& cur,uint64_t i,uint64_t o) {
            uint64_t offset;
            entry_t entry;
            _file.read(&entry,o+i*sizeof(entry_t),sizeof(entry_t));
            byteswap<SourceEndian>(&entry);
            offset=entry.count;
            if((entry.count*sizeof_type(entry.type)) <= sizeof(entry.data)) {
                offset=o+i*sizeof(entry_t)     // offset of current entry
                    +offsetof(entry_t,data);   // offset to data field
            } else {
                offset=entry.data;
            }
            cur.tags[(TiffTag)entry.tag]=data_t{(TiffType)entry.type,entry.count,offset};
        }

        template<typename T,Endian SourceEndian> vector<T> readtagdata_raw(const data_t& tag) {
            vector<T> buf(tag.nelem);
            _file.read(buf.data(), tag.offset, tag.nbytes());
            byteswapv<SourceEndian>(buf);
            return buf;
        }

        template<Endian SourceEndian,typename Tdst>
        vector<Tdst> readTagData(const data_t & tag) {
            vector<Tdst> v(tag.nelem);
            switch(tag.type) {
#define CASE(TagType,T) \
                case TagType: \
                    { auto data = readtagdata_raw<T,SourceEndian>(tag); \
                        copy(data.begin(), data.end(), v.begin()); \
                    } break;

                CASE(TIFFTYPE_BYTE,      uint8_t);
                CASE(TIFFTYPE_ASCII,     uint8_t);
                CASE(TIFFTYPE_SHORT,     uint16_t);
                CASE(TIFFTYPE_LONG,      uint32_t);
                CASE(TIFFTYPE_SBYTE,     int8_t);
                CASE(TIFFTYPE_SLONG,     int32_t);
                CASE(TIFFTYPE_FLOAT,     float);
                CASE(TIFFTYPE_DOUBLE,    double);
                CASE(TIFFTYPE_LONG8,     uint64_t);
                CASE(TIFFTYPE_SLONG8,    int64_t);
                CASE(TIFFTYPE_IFD8,      uint64_t);

                //Special cases
                case TIFFTYPE_RATIONAL: {
                    auto data = readtagdata_raw<uint32_t,SourceEndian>(tag);
                    copy(data.begin(), data.end(), v.begin());
                } break;
                case TIFFTYPE_SRATIONAL: {
                    auto data = readtagdata_raw<int32_t,SourceEndian>(tag);
                    copy(data.begin(), data.end(), v.begin());
                } break;
#undef CASE
                default:
                    throw runtime_error("Improper type for tag.  Or don't know how to handle.");
            }
            return v;
        }

        template<Endian SourceEndian,typename Tdst>
        Tdst readSingleValueFromTag(const data_t &tag) {
            auto v=readTagData(tag);
            return v.at(0);
		}

        template<Endian SourceEndian,typename Tdst>
        Tdst readSingleValueFromTag(const index_t &ifd,TiffTag tag_id) {
            auto it = ifd.tags.find(tag_id);
            if (it == ifd.tags.end()) {
                stringstream ss;
                ss << "Could not find tag " << tag_id;
                throw runtime_error(ss.str());
            }
            auto v=readTagData<SourceEndian,Tdst>(it->second);
            return v.at(0);
        }

        template<Endian SourceEndian,typename Tdst>
        Tdst readSingleValueFromTag(const index_t &ifd,TiffTag tag_id,const Tdst& default_value) {
            auto it = ifd.tags.find(tag_id);
            if (it == ifd.tags.end())
                return default_value;
            auto v=readTagData<SourceEndian,Tdst>(it->second);
            return v.at(0);
        }

        template<Endian SourceEndian>
        void readStripInfo(index_t & cur) {
            vector<uint64_t> offsets, nbytes;
            {
                auto it = cur.tags.find(TIFFTAG_STRIPBYTECOUNTS);
                if (it == cur.tags.end())
                    throw runtime_error("No strip offsets found.  This reader only reads Tiffs with strip data");
                nbytes = readTagData<SourceEndian,uint64_t>(it->second);
            }
            {
                auto it = cur.tags.find(TIFFTAG_STRIPOFFSETS);
                if (it == cur.tags.end())
                    throw runtime_error("No strip offsets found.  This reader only reads Tiffs with strip data");
                offsets = readTagData<SourceEndian,uint64_t>(it->second);
            }
            if (offsets.size() != nbytes.size())
                throw runtime_error("Number of StripOffsets and StripByteCounts disagree.");
            cur.data.resize(offsets.size());
            for (auto i = 0; i < cur.data.size();++i) {
                auto& d = cur.data[i];
                d.nbytes = nbytes[i];
                d.offset = offsets[i];
            }
        }

        template<Endian SourceEndian>
        void readBitsPerSample(index_t & cur) {
            cur.bits_per_sample = readSingleValueFromTag
                <SourceEndian,decltype(cur.bits_per_sample)>
                (cur, TIFFTAG_BITSPERSAMPLE);
        }

        template<Endian SourceEndian>
        nd_type readAndResolvePixelType(const index_t &ifd,unsigned bits_per_sample,unsigned samples_per_pixel) {
            TiffSampleFormat sample_format = (TiffSampleFormat)readSingleValueFromTag<SourceEndian,unsigned>(ifd, TIFFTAG_SAMPLEFORMAT,TIFFSAMPLEFORMAT_UNSIGNED);
            switch(sample_format) {
            case TIFFSAMPLEFORMAT_FLOATING:
                switch(bits_per_sample) {
                    case 32: return nd_f32;
                    case 64: return nd_f64;
                default:
                    goto Error;
                }
            case TIFFSAMPLEFORMAT_SIGNED:
                switch(bits_per_sample) {
                    case 8:  return nd_i8;
                    case 16: return nd_i16;
                    case 32: return nd_i32;
                    case 64: return nd_i64;
                default:
                    goto Error;
                }
            case TIFFSAMPLEFORMAT_UNSIGNED:
                switch(bits_per_sample) {
                    case 8:  return nd_u8;
                    case 16: return nd_u16;
                    case 32: return nd_u32;
                    case 64: return nd_u64;
                default:
                    goto Error;
                }
            case TIFFSAMPLEFORMAT_UNDEFINED:
            default:;
                goto Error;
            }
        Error:
            throw runtime_error("Could not resolve pixel type.");
		}

        /// Generic implementation used to build the index.
        ///
        /// Iterates over the ifds and build's an index that tells us where the
        /// interesting parts of the file are.
        ///
        template<typename nentries_t,typename offset_t,typename entry_t,Endian SourceEndian>
        void gen_index(){
            auto o=_first_ifd;
            nplanes = 0;
            // scan ifds for pointers
            while(o) {
                nentries_t nentries;
                _index.push_back(index_t());
                auto& cur=_index.back();
                _file.read(&nentries,o,sizeof(nentries));
                byteswap<SourceEndian>(&nentries);
                o+=sizeof(nentries);
                // scan ifd entries
                for(nentries_t i=0;i<nentries;++i)
                    readtag<entry_t,SourceEndian>(cur,i,o);
                readStripInfo<SourceEndian>(cur);
                readBitsPerSample<SourceEndian>(cur);
                if (nplanes == 0) { // special processing for first frame
                    w = readSingleValueFromTag<SourceEndian, decltype(w)>(cur, TIFFTAG_IMAGEWIDTH);
                    h = readSingleValueFromTag<SourceEndian, decltype(h)>(cur, TIFFTAG_IMAGELENGTH);
                    nchan = readSingleValueFromTag<SourceEndian,decltype(nchan)>(cur, TIFFTAG_SAMPLESPERPIXEL,1);
                    type = readAndResolvePixelType<SourceEndian>(cur,cur.bits_per_sample,nchan);
                }
                // For normal tiff this reads fewer bytes than sizeof(o), but
                // in that case the high bits of o are always 0 anyway.
                _file.read(&o,o+nentries*sizeof(entry_t),sizeof(offset_t));
                byteswap<SourceEndian,offset_t>((offset_t*)&o);
                ++nplanes;
                cur.next=o;
            }
        }
    };

} // vidrio::scanimage::tiff
} // vidrio::scanimage
} // vidrio


/*
TODO/FIXME

    TODO: Tag validation.
                Check internal consistency of tags.
                Check for (un)supported features.
                Check for valid tag type.

    TODO: refactoring.  names like index and strip and data etc.

    TODO: properly handle orientation

    TODO: TEST that endian-ness of image data is getting transformed
          properly.  It might just be getting done to a copy and not
          in-place.

    TODO: maybe separate out parts that assume the tiff is a stack (that
          all images are the same size).

          Also, validate assumption that all images in input
          are the same size

    TODO:  bits per sample is actually an array.  Properly we would unpack those
           different components bitwise into the output array.  Right now we act as
           if each sample has the same bitsPerSample (e.g. 8,8,8 for RGB) and that
           the number of bits aligns to bytes (is divisible by 8).

    TODO: add validation so we're at least up front about what we don't support

    TOOD: handle a "planar" PlanarConfiguration.  This mainly effects where data is
          copied when reading in the data.

    TODO: if we move to having roi's stored in different ifd's, then not all ifd's
          will be the same size.  The code will need to be adapted for that.

          Maybe the best approach will be to index ifd's based on (width,height).

    FIXME: remove data conversion warnings (from the tiff tag read single value function)

    FIXME: positively identify files with no metadata/metadata section
*/

/* NOTES

    1. The scanimage metadata/metadata section:

        At some point binary data was appended to the end of the tiff.  The size
        of the block was indicated by the last four bytes of the file.

        Unfortunately, there's no simple way to detect if a given Tiff is one of
        these files.  The method I'm using now is to try to detect the last byte
        referenced in the tiff file index.

        To compute the last byte, one properly needs to compare the last byte of
        the IFD section, and the last byte referenced by each tag.  Being lazy,
        I only look at the last byte of the image strip data.  For tiff's
        produced by ScanImage this should be fine, and I also imagine many other
        Tiff producers put their image data at the end, but it doesn't have to
        be this way.

        One day, I may have to do this the right way, but perhaps then I won't
        have to worry about the era wherein this particular version of the
        ScanImage Tiff file was written.

        We're altering the way we write the metadata section in the future so
        that it's can be positively identified.

	2. A reflection on my use of templates:

		I use templates for two things here.  (a) compile-time configuguration
		of the io dependency (file_t) and (b) generic programing to handle
		some of the different data-handling paths required when decoding tiffs.

		(b) was a good use of templates.

		(a) was bad.  I should just have declared the interface I wanted and
		then linked different implementations in based on the operating
		system.

		Templates have several downsides that I won't enumerate here, because
		they're obvious once you've used them a bit.  So they're worth avoiding.
		Link-time resolution of dependencies is what linking is for!  That
		would have turned out much cleaner.

*/
