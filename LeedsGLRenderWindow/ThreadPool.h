#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <queue>
#include <thread>
#include <mutex>
#include <functional>
#include <memory>
#include <pthread.h>
#include <condition_variable>

namespace TP {
    typedef std::function<void(void)> Task;
    typedef std::vector<std::function<void(void)> > TaskGroup;
    template<typename T>
    class TSafeQueue {
        // a thread safe queue
    private:
        std::mutex mut;
        std::queue<T> data_queue;

    public:
        TSafeQueue() {
        }

        TSafeQueue(TSafeQueue const &other) {
            std::lock_guard<std::mutex> lock(other.mut);
            data_queue = other.data_queue;
        }

        void push(T value) {
            std::lock_guard<std::mutex> lock(mut);
            data_queue.push(value);
        }

        std::shared_ptr<T> pop() {
            std::lock_guard<std::mutex> lock(mut);
            std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
            data_queue.pop();
            return res;
        }

        std::shared_ptr<T> try_pop() {
            std::lock_guard<std::mutex> lock(mut);
            if (data_queue.empty())
                return std::shared_ptr<T>();
            std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
            data_queue.pop();
            return res;
        }

        bool size() {
            std::lock_guard<std::mutex> lock(mut);
            return data_queue.size();
        }

        bool empty() {
            std::lock_guard<std::mutex> lock(mut);
            return data_queue.empty();
        }

        void clear() {
            std::lock_guard<std::mutex> lock(mut);
            while (!data_queue.empty())data_queue.pop();
        }
    };

    class WorkerThread {
        // a worker thread
    private:
        TSafeQueue<Task> queue;
        std::thread t;
        std::mutex mux;
        std::condition_variable thread_returned;
        std::thread::id tid;
        bool is_running;
        bool thread_stopped;

        bool isRunning() {
            std::unique_lock<std::mutex> lock(mux);
            return is_running;
        }

        void isRunning(bool running) {
            std::lock_guard<std::mutex> lock(mux);
            is_running = running;
            if (running) {
                thread_stopped = false;
            }
        }

    public:
        WorkerThread() {
            is_running = false;
            thread_stopped = true;
        };

        ~WorkerThread() {
            queue.clear();
            stop();
        }

        void run() {
            isRunning(true);
            t = std::move(std::thread(&WorkerThread::mainLoop, this));
            tid = t.get_id();
            t.detach();
        }

        void stop() {
            isRunning(false);
            std::unique_lock<std::mutex> lock(mux);
            // wait for the thread finish all works
            thread_returned.wait(lock, [this] {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                return thread_stopped;
            });
        }

        void insertTask(Task &task) {
            queue.push(task);
        }

        int remainTasks() {
            return queue.size();
        }

        // main loop
        void mainLoop() {
            do {
                if (isRunning() && !queue.empty()) {
                    Task task = *(queue.pop());
                    task();
                } else {
                    // else sleep for a while
                    std::this_thread::sleep_for(std::chrono::milliseconds (1));
                }
            } while (isRunning());
            std::unique_lock<std::mutex> lk(mux);
            thread_stopped = true;
            thread_returned.notify_one();
        }
    };
    class ThreadPool {
    private:
        TSafeQueue<Task> queue;
        std::vector<WorkerThread *> workers;
        std::condition_variable is_empty;
        unsigned int current_index;
        unsigned int max_index;

    public:
        explicit ThreadPool(unsigned int count) {
            max_index = count;
            current_index = 0;
            for (unsigned int i = 0; i < count; i++) {
                auto *worker = new WorkerThread();
                worker->run();
                workers.push_back(worker);
            }
        };

        ~ThreadPool() {
            for (auto w: workers) {
                w->stop();
                delete w;
            }
        }

        void doAsync(Task task) {
            current_index = (current_index + 1) % max_index;
            workers[current_index]->insertTask(task);
        }

        void doSync(Task task) {
            std::mutex mux;
            mux.lock();
            doAsync([&mux, &task]() {
                if (task) {
                    task();
                }
                mux.unlock();
            });
            mux.lock();
        }

        void syncGroup(TaskGroup &tasks, int batch_size = -1) {
            if (batch_size == -1) {
                // default size, alloc tasks with avg count to every worker
                batch_size = int(tasks.size() / workers.size());
            }
            batch_size = std::max(batch_size, 0);
            if (batch_size == 0) {
                // execute directly
                for (auto &task : tasks) {
                    task();
                }
                return;
            }
            TaskGroup real_tasks;
            if (batch_size >= 1) {
                // alloc tasks to all workers
                int real_size = int(tasks.size()) / batch_size + 1;
                real_tasks.resize(real_size);
                for (unsigned int i = 0; i < real_size; i++) {
                    real_tasks[i] = [&tasks, i, batch_size]() {
                        for (int j = 0; j < batch_size; j++) {
                            int index = i * batch_size + j;
                            if (index < tasks.size()) {
                                tasks[index]();
                            }
                        }
                    };
                }
            }
            // use count to sync all tasks
            int count = real_tasks.size();
            std::mutex mx;
            for (auto task: real_tasks) {
                doAsync([task, this, &mx, &count]() {
                    task();
                    std::unique_lock<std::mutex> lock(mx);
                    if ((--count) == 0) {
                        is_empty.notify_one();
                    }
                });
            }
            std::unique_lock<std::mutex> lock(mx);
            is_empty.wait(lock, [&count] {
                return count == 0;
            });
        }
    }; // ThreadPool
}; // TP

#endif // THREAD_POOL_H