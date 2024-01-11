#pragma once

#define SPLASH_VER "1.9.0"

inline void SPLASH_VER_PRINT(std::ostream& oss) {
	oss << "Version: " << SPLASH_VER << "\n";
}
