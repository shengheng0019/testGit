#pragma once

#include <mutex>
#include <queue>

// Thread safe implementation of a Queue using an std::queue
template<class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class TaskQueue {
private:
        std::priority_queue<T, Container, Compare> m_queue;
        std::mutex m_mutex;
public:
        TaskQueue() {

        }

        TaskQueue(TaskQueue& other) {
                //TODO:
        }

        ~TaskQueue() {

        }


        bool empty() {
                std::unique_lock<std::mutex> lock(m_mutex);
                return m_queue.empty();
        }

        int size() {
                std::unique_lock<std::mutex> lock(m_mutex);
                return m_queue.size();
        }

        void enqueue(T& t) {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_queue.push(t);
        }

        bool dequeue(T& t) {
                std::unique_lock<std::mutex> lock(m_mutex);

                if (m_queue.empty()) {
                        return false;
                }
                t = std::move(m_queue.top());

                m_queue.pop();
                return true;
        }
};