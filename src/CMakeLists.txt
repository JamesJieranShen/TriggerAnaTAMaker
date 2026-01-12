file(GLOB SRC_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
)

add_executable(MakeTriggerActivity ${SRC_FILES})

target_link_libraries(MakeTriggerActivity PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::TreePlayer ROOT::Physics)
target_link_libraries(MakeTriggerActivity PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json)

target_include_directories(MakeTriggerActivity PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

install(TARGETS MakeTriggerActivity)

