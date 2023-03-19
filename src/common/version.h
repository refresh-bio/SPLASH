#pragma once

#define NOMAD_VER "2.0.0"

inline void NOMAD_VER_PRINT(std::ostream& oss) {
	oss << "Version: " << NOMAD_VER << "\n";
}
