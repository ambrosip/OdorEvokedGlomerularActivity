#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <string>

using namespace std;

string ProgramName() {
    char buf[1024]={0};
    // get the fully qualified path name to the program
    DWORD n=GetModuleFileName(GetModuleHandle(0),buf,sizeof(buf));
    // get rid of the path and the extension
    unsigned ext=n;
    while(ext>0 && buf[--ext]!='.');
    if(ext==0) ext=n;
    buf[ext]='\0';
    unsigned pathsep=ext;
    while(pathsep>0 && buf[--pathsep]!='\\');
    return string(buf+pathsep+1,ext-pathsep-1);
}

