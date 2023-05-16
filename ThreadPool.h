#pragma once

#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>
#include <iostream>
#include <utility>
#include <functional>
#include <chrono>
#include <thread>

#include "TaskQueue.h"

typedef std::pair<int,std::function<void()>> QuenuEvent;
struct CustomComparator {
  bool operator()(const QuenuEvent& left, const QuenuEvent& right) const {
    return left.first > right.first;
  }
};

class ThreadPool
{
private:
        class ThreadWorker
        {
        private:
                int m_id;
                ThreadPool* m_pool;
        public:
                ThreadWorker(ThreadPool* pool, const int id)
                        : m_pool(pool), m_id(id) {
                }

                void operator()() {
                        QuenuEvent event;
                        bool dequeued;
                        while (!m_pool->m_shutdown) {
                                {
                                        std::unique_lock<std::mutex> lock(m_pool->m_conditional_mutex);
                                        if (m_pool->m_queue.empty()) {
                                                m_pool->m_conditional_lock.wait(lock);
                                        }
                                        dequeued = m_pool->m_queue.dequeue(event);
                                }
                                if (dequeued) {
                                        event.second();
                                        std::this_thread::sleep_for(std::chrono::milliseconds(200));
                                }
                        }
                }
        };

        bool m_shutdown;

        TaskQueue<QuenuEvent , std::vector<QuenuEvent>, CustomComparator> m_queue;
        std::vector<std::thread> m_threads;
        std::mutex m_conditional_mutex;
        std::condition_variable m_conditional_lock;

public:
        ThreadPool(const int n_threads)
                : m_threads(std::vector<std::thread>(n_threads)), m_shutdown(false) {
        }

        ThreadPool(const ThreadPool&) = delete;
        ThreadPool(ThreadPool&&) = delete;

        ThreadPool& operator=(const ThreadPool&) = delete;
        ThreadPool& operator=(ThreadPool&&) = delete;

        // Inits thread pool
        void init() {
                for (int i = 0; i < m_threads.size(); ++i) {
                        m_threads[i] = std::thread(ThreadWorker(this, i));
                }
        }

        // Waits until threads finish their current task and shutdowns the pool
        void shutdown() {
                m_shutdown = true;
                m_conditional_lock.notify_all();

                for (int i = 0; i < m_threads.size(); ++i) {
                        if (m_threads[i].joinable()) {
                                m_threads[i].join();
                        }
                }
        }

        // Submit a function to be executed asynchronously by the pool
        template<typename T, typename F, typename...Args>
        auto submit(const T t, F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
                // Create a function with bounded parameters ready to execute
                std::function<decltype(f(args...))()> func = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
                // Encapsulate it into a shared ptr in order to be able to copy construct / assign 
                auto task_ptr = std::make_shared<std::packaged_task<decltype(f(args...))()>>(func);

                // Wrap packaged task into void function
                std::function<void()> wrapper_func = [task_ptr]() {
                        (*task_ptr)();
                };

                int p = static_cast<int>(t);

                QuenuEvent event = std::make_pair(p, wrapper_func);

                // Enqueue generic wrapper function
                m_queue.enqueue(event);

                // Wake up one thread if its waiting
                m_conditional_lock.notify_one();

                // Return future from promise
                return task_ptr->get_future();
        }
};
