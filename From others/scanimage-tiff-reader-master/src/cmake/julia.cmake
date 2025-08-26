install(
    FILES julia/ScanImageTiffReader.jl
    DESTINATION share/julia
    COMPONENT julia)


install(
    TARGETS ScanImageTiffReaderAPI
    DESTINATION share/julia
    COMPONENT julia
)
