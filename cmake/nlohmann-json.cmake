if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

set(NLOHMANNJSON_VERSION "v3.9.1")

# include(FetchContent)
# FetchContent_Declare(
#     nlohmann_json
#     URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
#     URL_HASH SHA256=87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014
# )

include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
    URL_HASH SHA256=6bea5877b1541d353bd77bdfbdb2696333ae5ed8f9e8cc22df657192218cad91
)
FetchContent_MakeAvailable(nlohmann_json)

add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE
    ${nlohmann_json_SOURCE_DIR}/include)
add_library(nlohmann::json ALIAS nlohmann_json)
