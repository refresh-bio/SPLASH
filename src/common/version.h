#pragma once

#define SPLASH_VER "2.0.3"

inline void SPLASH_VER_PRINT(std::ostream& oss) {
	oss << "Version: " << SPLASH_VER << "\n";
}
