#pragma once

#define SPLASH_VER "2.11.2"

inline void SPLASH_VER_PRINT(std::ostream& oss) {
	oss << "splash version: " << SPLASH_VER << "\n";
}
