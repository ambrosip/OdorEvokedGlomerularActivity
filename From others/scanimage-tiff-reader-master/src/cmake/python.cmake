install(FILES python/ScanImageTiffReader.py
        DESTINATION share/python
        COMPONENT python)
install(
    TARGETS ScanImageTiffReaderAPI
    DESTINATION share/python
    COMPONENT python
)
