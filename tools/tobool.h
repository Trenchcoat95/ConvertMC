#ifndef tobool_h
#define tobool_h

#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cctype>

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

#endif