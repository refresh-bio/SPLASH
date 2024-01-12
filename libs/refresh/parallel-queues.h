#ifndef _REFRESH_PARALLEL_QUEUES
#define _REFRESH_PARALLEL_QUEUES

#include "parallel-queues-common.h"
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <map>

namespace refresh {

	const uint32_t REFRESH_BUILD_PARALLEL_QUEUES = 1;

	template<typename T>
	class circular_queue
	{
		std::vector<T> data;
		bool full = false;
		size_t start = 0;
		size_t end = 0;
	public:
		circular_queue(size_t size) :
			data(size)
		{

		}
		bool empty() const
		{
			return start == end && !full;
		}
		void emplace(T&& elem)
		{
			data[end] = std::move(elem);
			end = (end + 1ul) % data.size();
			full = end == start;
		}

		bool is_full() const
		{
			return full;
		}

		T& front()
		{
			return data[start];
		}

		void pop()
		{
			start = (start + 1ul) % data.size();
			full = false;
		}
	};

	template<typename T>
	class limited_stl_queue : public std::queue<T>
	{
		size_t max_size;
	public:
		limited_stl_queue(size_t size) :
			max_size(size)
		{

		}
		bool is_full() const
		{
			return this->size() >= max_size;
		}
	};

	//implements thread safe queue
	template<typename T, typename QUEUE_T = circular_queue<T>>
	class parallel_queue
	{
		QUEUE_T q;
		bool is_completed = false;
		bool canceled = false;
		size_t n_writers;
		std::string name;
		std::mutex mtx;
		std::condition_variable cv_push;
		std::condition_variable cv_pop;

		IQueueObserver* queue_observer;

		void push_impl(T&& elem)
		{
			bool was_empty = q.empty();
			q.emplace(std::move(elem));

			if (was_empty)
				cv_pop.notify_one();

			if (queue_observer)
				queue_observer->notify_pushed();
		}

	public:
		parallel_queue(
			size_t size,
			size_t n_writers = 1,
			const std::string& name = "",
			IQueueObserver* queue_observer = nullptr)
			:
			q(size),
			n_writers(n_writers),
			name(name),
			queue_observer(queue_observer)
		{
#ifdef REFRESH_PROFILE_QUEUES
			if (name.empty())
				this->name = refresh_profile_queues_create_queue_name();
			if (!queue_observer)
				this->queue_observer = refresh_profile_queues_create_observer();
#endif
			if (this->queue_observer)
				this->queue_observer->set_queue_params(this->name, size);
		}

		//returns false if pushing thread should stop
		bool push_or_cancel(T&& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();
			std::unique_lock<std::mutex> lck(mtx);

			cv_push.wait(lck, [this] {return !q.is_full() || canceled; });

			if (queue_observer)
				queue_observer->notify_wait_on_push_time(std::chrono::high_resolution_clock::now() - time_start);

			if (canceled)
				return false;

			push_impl(std::move(elem));
			return true;
		}

		void push(T&& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();
			std::unique_lock<std::mutex> lck(mtx);

			cv_push.wait(lck, [this] {return !q.is_full(); });

			if (queue_observer)
				queue_observer->notify_wait_on_push_time(std::chrono::high_resolution_clock::now() - time_start);

			push_impl(std::move(elem));
		}

		bool pop(T& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();
			std::unique_lock<std::mutex> lck(mtx);

			cv_pop.wait(lck, [this] {
				return !q.empty() || is_completed;
				});
			if (queue_observer)
				queue_observer->notify_wait_on_pop_time(std::chrono::high_resolution_clock::now() - time_start);

			if (is_completed && q.empty())
				return false;

			bool was_full = q.is_full();

			elem = std::move(q.front());
			q.pop();

			if (was_full)
				cv_push.notify_one();

			if (queue_observer)
				queue_observer->notify_popped();
			return true;
		}

		void mark_completed()
		{
			std::lock_guard<std::mutex> lck(mtx);
			--n_writers;
			if (!n_writers)
			{
				is_completed = true;
				cv_pop.notify_all();
			}
		}

		void cancel()
		{
			std::lock_guard<std::mutex> lck(mtx);
			canceled = true;
			cv_push.notify_all();
		}
	};

	//implements thread safe priority queue
	template<typename T>
	class parallel_priority_queue
	{
		std::map<uint64_t, T> map_data;

		size_t max_size;
		bool is_completed = false;
		size_t n_writers;

		std::string name;
		IQueueObserver* queue_observer;

		std::mutex mtx;
		std::condition_variable cv_push;
		std::condition_variable cv_pop;
		uint64_t current_priority{};

	public:
		parallel_priority_queue(
			size_t size,
			uint64_t n_writers = 1,
			const std::string& name = "",
			IQueueObserver* queue_observer = nullptr) :
			max_size(size),
			n_writers(n_writers),
			name(name),
			queue_observer(queue_observer)
		{
#ifdef REFRESH_PROFILE_QUEUES
			if (name.empty())
				this->name = refresh_profile_queues_create_queue_name();
			if (!queue_observer)
				this->queue_observer = refresh_profile_queues_create_observer();
#endif
			if (this->queue_observer)
				this->queue_observer->set_queue_params(this->name, size);

		}
		void push(uint64_t priority, T&& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();

			std::unique_lock<std::mutex> lck(mtx);
			cv_push.wait(lck, [this, &priority] {
				return map_data.size() < max_size - 1 ||
					priority == current_priority;
				});

			if (queue_observer)
				queue_observer->notify_wait_on_push_time(std::chrono::high_resolution_clock::now() - time_start);

			map_data.emplace(priority, std::move(elem));

			if (priority == current_priority)
				cv_pop.notify_all();

			cv_push.notify_all();

			if (queue_observer)
				queue_observer->notify_pushed();
		}

		bool pop(T& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();

			std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this] {
				return is_completed || (!map_data.empty() && map_data.begin()->first == current_priority);
				});

			if (queue_observer)
				queue_observer->notify_wait_on_pop_time(std::chrono::high_resolution_clock::now() - time_start);

			if (is_completed && map_data.empty())
				return false;

			elem = std::move(map_data.begin()->second);
			map_data.erase(map_data.begin());

			++current_priority;

			cv_push.notify_all();

			if (queue_observer)
				queue_observer->notify_popped();

			return true;
		}

		//mkokot_TODO: implement Cancel and PushOrCancel

		void mark_completed()
		{
			std::lock_guard<std::mutex> lck(mtx);
			--n_writers;
			if (!n_writers)
			{
				is_completed = true;
				cv_pop.notify_all();
			}
		}
	};


	//implements thread safe queue
	/*
	It is assumed that there is exactly one producer who can steal idle threads
	Helps to implement work-stealing-like approach
	*/
	template<typename T, typename QUEUE_T = circular_queue<T>>
	class parallel_queue_pop_waiting
	{
		QUEUE_T q;
		bool is_completed = false;
		std::mutex mtx;
		std::condition_variable cv_push;
		std::condition_variable cv_pop;
		uint32_t n_pop_waiting{};

		IQueueObserver* queue_observer;
		std::string name;

	public:
		parallel_queue_pop_waiting(
			size_t size,
			const std::string& name = "",
			IQueueObserver* queue_observer = nullptr) :
			q(size),
			name(name),
			queue_observer(queue_observer)
		{
#ifdef REFRESH_PROFILE_QUEUES
			if (name.empty())
				this->name = refresh_profile_queues_create_queue_name();
			if (!queue_observer)
				this->queue_observer = refresh_profile_queues_create_observer();
#endif
			if (this->queue_observer)
				this->queue_observer->set_queue_params(this->name, size);

		}
		size_t get_n_waiting_on_pop() const
		{
			return n_pop_waiting;
		}
		void push(T&& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();

			std::unique_lock<std::mutex> lck(mtx);

			cv_push.wait(lck, [this] {return !q.is_full(); });

			if (queue_observer)
				queue_observer->notify_wait_on_push_time(std::chrono::high_resolution_clock::now() - time_start);

			bool was_empty = q.empty();
			q.emplace(std::move(elem));

			if (was_empty)
				cv_pop.notify_one();

			if (queue_observer)
				queue_observer->notify_pushed();
		}

		bool pop(T& elem)
		{
			auto time_start = std::chrono::high_resolution_clock::now();

			std::unique_lock<std::mutex> lck(mtx);
			bool first = true;
			cv_pop.wait(lck, [this, &first] {
				if (first)
					++n_pop_waiting;
				first = false;
				return !q.empty() || is_completed;
				});

			if (queue_observer)
				queue_observer->notify_wait_on_pop_time(std::chrono::high_resolution_clock::now() - time_start);

			--n_pop_waiting;

			if (is_completed && q.empty())
				return false;

			bool was_full = q.is_full();
			elem = std::move(q.front());
			q.pop();

			if (was_full)
				cv_push.notify_one();

			if (queue_observer)
				queue_observer->notify_popped();

			return true;
		}

		void mark_completed()
		{
			std::lock_guard<std::mutex> lck(mtx);
			is_completed = true;
			cv_pop.notify_all();
		}
	};

} //namespace refresh
#endif // _REFRESH_PARALLEL_QUEUES
