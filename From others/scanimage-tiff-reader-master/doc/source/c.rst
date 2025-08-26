C Library
=========

The C API is provided as a static or shared library. It requires no
dependencies.  We try to statically link wherever possible.

Open/Close
~~~~~~~~~~

Use :api:`ScanImageTiffReader_Open` to open a file.  When
you're done call :api:`ScanImageTiffReader_Close` to close any open
file handles and cleanup any other resources in use.

Error handling
~~~~~~~~~~~~~~

When errors occur, functions usually return a ``0``.  Any error messages will be
contained in the ``log`` field of the :api:`ScanImageTiffReader` context object.
Normally the ``log`` string is ``NULL``; only after an error does it point to
a null-terminated string.

Image data
~~~~~~~~~~

To read the image data,

    1. Determine the size of the buffer you'll need to hold the raw data using
       :api:`ScanImageTiffReader_GetDataSizeBytes`.
    2. Allocate the memory you need.
    3. Read using :api:`ScanImageTiffReader_GetData`.  This loads the entire
       volume into memory.  There's no option for reading just a subset of the
       data.

One nice thing about :api:`ScanImageTiffReader_GetData` is that many asynchronous
read requests are made at once.  This typically means that read commands will
be coallesced and queued by the operating system so that very high read
bandwidths can be acheived.

Metadata
~~~~~~~~

Reading metadata follows a similar pattern.  To read the frame-invariant metadata:

    1. Determine the size of the buffer you need using
       :api:`ScanImageTiffReader_GetMetadataSizeBytes`.
    2. Read it using :api:`ScanImageTiffReader_GetMetadata`.

The frame-varying metadata is read using a similar pattern.  There are two
sets of functions for reading the frame varying data:

    * :api:`ScanImageTiffReader_GetImageDescription` reads the image description
      for a single frame.
    * :api:`ScanImageTiffReader_GetAllImageDescriptions` reads the image
      descriptions for all frames into a single buffer.


