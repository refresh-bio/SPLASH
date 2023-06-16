#pragma once

#define NOMAD_VER "2.1.4"

inline void NOMAD_VER_PRINT(std::ostream& oss) {
	oss << "Version: " << NOMAD_VER << "\n";
}
