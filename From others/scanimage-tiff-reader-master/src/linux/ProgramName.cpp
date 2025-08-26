#include <string>
#include <exception>
#include <stdexcept>
#include <cstring>
#include <unistd.h>

using namespace std;

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

string ProgramName() {
    char buf[1024]={0};
#ifdef __APPLE__
    {
        unsigned sz=sizeof(buf);
        if(_NSGetExecutablePath(buf,&sz)!=0)
          throw runtime_error("Failed to get name of executable from _NSGetExecutablePath.");
        return string(buf);
    }
#else

    ssize_t len=readlink("/proc/self/exe",buf,sizeof(buf)-1);
    if(len!=-1) {
        buf[len]='\0';
        return string(buf);
    }
    throw runtime_error("Failed to get name of executable from /proc/self/exe.");
#endif

}
