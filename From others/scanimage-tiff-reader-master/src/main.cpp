#include "tiff.reader.h"
#include "platform.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <map> // used for argument processing

using namespace std;
using namespace vidrio;

static const char* AUTHOR="Vidrio Technologies <support@vidriotech.com>";

/// \returns the name of the current executable.
extern string ProgramName();

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)

static void headline() {
	string name=ProgramName();
	printf("%s %s-%s by %s\n",name.c_str(),STR(GIT_TAG),STR(GIT_HASH),AUTHOR);
}

/// Configure io implementation to use.
/// \see vidrio::platform::io
#if defined(_MSC_VER) || defined(WIN32)
    // about a factor of 2 faster to use native file io rather than
    // the standard library (Release build, Win 10 x64, March 2016)
    #include "win32/io.h"
    typedef platform::io<platform::win32::io_impl> io;
#else
    #include "linux/io.h"
    //typedef platform::io<platform::standard::io_impl> io;
    typedef platform::io<platform::linux_pread::io_impl> io;
#endif

static int allheaders(int argc, char**argv) {
    scanimage::tiff::reader<io> f(argv[0]);
    ostream *out = &cout;
    ofstream fout;
    if(argc==2) {
        fout.open(argv[1]);
        out = &fout;
    }
    
    unsigned frame=0;
    for(auto h:f.headers())
        (*out) << "[FRAME " << frame++ << " HEADER]" << endl 
               << h << endl << endl;
    return argc; // returns number of consumed arguments
}

static int oneheader(int argc, char**argv) {
    
    auto i=stoi(argv[0]);
    scanimage::tiff::reader<io> f(argv[1]);
    ostream *out = &cout;
    ofstream fout;
    if(argc==3) {
        fout.open(argv[2]);
        out = &fout;
    }    
    try {
        (*out) << f.headers().at(i) << endl;
    } catch(out_of_range e) {
        // FIXME: Say some subset of frames lack image descriptors for 
        //        whatever reason.  Then the correspondance between frames
        //        and headers is weird.  Might prefer to just return nothing.
                 
        ;// do nothing
    }
    return 3; // returns number of consumed arguments
}

static int metadata(int argc, char**argv) {
    scanimage::tiff::reader<io> f(argv[0]);
    ostream *out = &cout;
    ofstream fout;
    if(argc==2) {
        fout.open(argv[1]);
        out = &fout;
    }
    (*out) << f.metadata() << endl;
    return 2; // returns number of consumed arguments
}


static int help(int argc, char**argv);

static int image(int argc, char**argv) {
    help(1, argv-1);
    return argc;
}

#include<cmath>
#define countof(e) (sizeof(e)/sizeof(*(e)))
string HumanReadable(uint64_t t) {
    double decades = log((double)t)/log(2.0);
    const char* designator[] = { " ", " kB", " MB", " GB", " TB", " PB" };
    int i = (int) floor(decades / 10);
    i = (i < 0) ? 0 : (i>countof(designator)) ? countof(designator) : i;
    stringstream s;
    s << ((double)t / (double) (1ULL << (10 * i))) << designator[i];
    return s.str();
}

static int image_nbytes(int argc, char**argv) {
    scanimage::tiff::reader<io> f(argv[0]);
    cout << HumanReadable(f.bytesof_data()) << endl;
    return argc;
}

class scoped_timer {
    decltype(chrono::high_resolution_clock::now()) t0_;
    string label_;
public:
    scoped_timer(string label) : label_(label) {
        t0_=chrono::high_resolution_clock::now();
    }
    void now() const {
        auto t1 = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> t_ms = t1 - t0_;
        cerr << label_ << " Elapsed: " << t_ms.count() << "ms" << endl;
    }
};

static int rawvol(int argc, char**argv) {
    scoped_timer t0("Open");
    scanimage::tiff::reader<io> f(argv[0]);
    t0.now();
    scoped_timer t1("Read");
    auto v = f.data();
    t1.now();
    {
        scoped_timer t1("Write");
        ofstream out(argv[1], ios::binary);
        out.write(v.data(), v.size());
        t1.now();
    }
    return argc;
}

static int shape(int argc, char**argv) {
    scanimage::tiff::reader<io> f(argv[0]);
    const char *typenames[] = { "u8","u16","u32","u64","i8","i16","i32","i64","f32","f64" };
    nd shape = f.shape();
    cout << "Shape: " << shape.dims[0];
    for (auto i = 1U; i < shape.ndim;++i) {
        cout << " x " << shape.dims[i];
    }
    cout << " @ " << typenames[shape.type] << endl;
    return argc;
}

/// \brief Data for defining command line arguments and options.
/// \see Commands
/// 
/// Notice it's got a sort of recursive structure.  That's not really used
/// at the moment.
struct command {
    list<string> args,
                 optargs;                           ///< required fields that follow the command
    string description;                             ///< Help string describing what the command does
    function<int(int argc, char** argv)> handler;   ///< Called when the command is found with no options. \returns number  of handled arguments
    map<string,command> options;                    ///< Optional fields modifying the command 
};

static map<string, command> Commands{
    {"descriptions",{
        {"input file"},{"output file"},
         "Extract the contents of a the image description tag(s). "
         "If you don't specify which frame, then all the image descriptions "
         "are output.",
          &allheaders,
          {{"--frame",{
              {"int"},{},
              "Extract the description for just the specified frame.",
              &oneheader
          }}}
    }},
    {"image",{
        {},{},
        "Do stuff with image data.",
        &image,
        {{"bytes",{
            {"input file"},{},
            "Return the byte size of the raw image data.",
            &image_nbytes
        }},{"raw",{
            {"input file","output file"},{},
            "Save raw bytes to a file for the entire volume.",
            &rawvol
        }},{"shape",{
            {"input file"},{},
            "Print the shape and type of the volume in the tiff stack.",
            &shape
        }}
        }
    }},
    {"metadata",{
        {"input file"},{"output file"},
         "Extract the ScanImage metadata section from the file.",
         &metadata
    }},
    {"help",{
        {},{"sub-command"},
        "Print detailed help for the specified sub-command.",
        &help
    }}
};

#include <iomanip>

static void usage_cmd(const string &name,const command &cmd) {
    cerr << '\t' << setw(15) << name <<'\t'<<cmd.description << endl;
}

static void usage_detailed_help(const string &context, const string &name,const command &cmd) {
    cerr << context << ' ' << name << ' ';
    if (cmd.options.size())
        cerr << "[options] ";
    for (auto &arg : cmd.args)
        cerr << '<' << arg << "> ";
    for (auto &arg : cmd.optargs)
        cerr << '[' << arg << "] ";
    cerr << endl;

    cerr << '\t' << cmd.description << endl << endl;
    
    for (auto &opt : cmd.options) {
        cerr << '\t' << opt.first << " ";
        for (auto &arg : opt.second.args)
            cerr << '<' << arg << "> ";
        cerr << '\t' << opt.second.description << endl << endl;
    }
}

static void usage() {	
    cerr << "Usage: " << ProgramName() << " <command> [<args>]" << endl
         << endl 
         << "These are the available commands: " << endl;
    for (auto &e : Commands)
        usage_cmd(e.first,e.second);
    cerr << endl;
}

static int help(int argc, char**argv)
{
    const char* cmd=argc?argv[0]:"help";
    auto el=Commands.find(cmd);
    if (el != Commands.end()) {
        usage_detailed_help(ProgramName(), el->first, el->second);
        if (el->first == "help") {
            cout << "\tAvailable sub-commands are:" << endl;
            for (auto &cmd : Commands)
                cout << "\t\t" << cmd.first << endl;
        }
    }
    else
        throw invalid_argument("Cold not find help for the specified command");
    return  1; // returns number of consumed arguments
}

static void parseargs(int argc, char**argv) {
    // skip the first one: it's something like the command used to call this program
    argc--;
    argv++;
    bool any_cmd = false;
    while (argc) {
        auto pcmd = Commands.find(*argv);
        if (pcmd != Commands.end()) {
            any_cmd = true;
            auto cmd = pcmd->second;
            argc--; argv++; // discard the command name
            // check for options
            auto any_opt = false;
            for (auto i = 0; i < argc; ++i) {
                auto popt = cmd.options.find(*argv);
                if(popt!=cmd.options.end()) {
                    auto opt = popt->second;
                    any_opt = true;
                    argc--; argv++;         // discard the option flag                    
                    if (argc < opt.args.size()) {
                        stringstream ss;
                        ss << "Expected arguments to " << pcmd->first << " " << popt->first;
                        throw invalid_argument(ss.str());
                    }
                    auto n=opt.handler(argc, argv);// got a match.  must modify argc,argv to consume option
                    argc -= (int) n;
                    argv += (int) n;
                }
            }
            if (!any_opt) {
                // no options found - call default_callback
                if (argc < cmd.args.size() || argc > cmd.args.size()+cmd.optargs.size()) {
                    help(1, argv-1);
                    stringstream ss;
                    ss << "Expected arguments to " << pcmd->first;
                    throw runtime_error(ss.str());
                }
                auto n=cmd.handler(argc, argv);// must modify argc,argv to consume option
                argc -= (int) n;
                argv += (int) n;
            }
        }
        break;
    }
    if (!any_cmd) {
        throw invalid_argument("Did not recognize command.");
    }
}

#include <chrono>
#define NREPEATS (1)

int main(int argc,char* argv[]){
    try{
		headline();
        for (int i = 0; i < NREPEATS; ++i) {
            auto t0 = chrono::high_resolution_clock::now();

            parseargs(argc, argv);            

            auto t1 = chrono::high_resolution_clock::now();
            chrono::duration<double, milli> t_ms = t1 - t0;
            cerr << "Elapsed: " << t_ms.count() << "ms" << endl;
        }
        return 0;
    } catch(invalid_argument &e) {
        usage();
        cerr << e.what() << endl;
        return 1;
    } catch(exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 1;
}

/// \file main.cpp
/// Entry function for command line utility.
/// Extracts the image descriptors and metadata from a ScanImage tiff.


// TODO: make this so the options don't have to come before the arguments.
//       this requires making argv into a list and then manipulating that,
//       which is probably safer anyway
