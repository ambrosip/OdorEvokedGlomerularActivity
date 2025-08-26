
#include <tiff.reader.api.h>

int main(int argc,char*argv[]) {
    ScanImageTiffReader reader=ScanImageTiffReader_Open("E:/Downloads/example.tif");
    
    printf("Number of frames: %d\n", (int) ScanImageTiffReader_GetFrameCount(&reader));
    
    size_t nbytes=ScanImageTiffReader_GetFrameIntervalSizeBytes(&reader,100,200);
    void *buf=malloc(nbytes);
    if(!buf) {
        printf("Could not allocate %llu bytes\n",nbytes);
        exit(1);
    }

    ScanImageTiffReader_GetFrameInterval(&reader,100,200,buf,nbytes);

    if(reader.log) {
        printf("Error\n\t%s\n",reader.log);
        exit(1);
    }

    ScanImageTiffReader_Close(&reader);
    return 0;
}