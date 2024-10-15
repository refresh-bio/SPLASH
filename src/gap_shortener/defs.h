#pragma once

#include <refresh/parallel_queues/lib/parallel-queues.h>
#include "unsafe_vector.h"
#include <cinttypes>
#include <memory>

constexpr uint32_t reads_pack_size = 2 << 18;

using read_t = std::string;
using read_pack_t = std::vector<read_t>;

constexpr uint32_t pseudoreads_pack_size = 2 << 25;

using pseudoreads_pack_t = unsafe_vector<uint8_t>;