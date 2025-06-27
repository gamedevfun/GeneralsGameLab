set(GLTF2W3D_SRC
    "gltf2w3d.cpp"
)

add_executable(core_gltf2w3d WIN32)
set_target_properties(core_gltf2w3d PROPERTIES OUTPUT_NAME gltf2w3d)

target_sources(core_gltf2w3d PRIVATE ${GLTF2W3D_SRC})

target_include_directories(core_gltf2w3d PRIVATE
    "."
    "../WW3D/max2w3d"
    "../WW3D/pluglib"
    "../../Libraries/Source/WWVegas"
    "../../Libraries/Source/WWVegas/WWLib"
)

target_link_libraries(core_gltf2w3d PRIVATE
    core_config
    corei_always
)

if(WIN32 OR "${CMAKE_SYSTEM}" MATCHES "Windows")
    target_link_options(core_gltf2w3d PRIVATE /subsystem:console)
endif()
