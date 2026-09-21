#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <memory>

// ---------- 线程安全任务队列 ----------
template<typename T>
class TaskQueue {
public:
    void push(const T& task) {
        std::lock_guard<std::mutex> lock(mtx_);
        queue_.push(task);
        cv_.notify_one();
    }

    bool pop(T& task) {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [this]{ return !queue_.empty() || finished_; });
        if (queue_.empty()) return false;
        task = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void finish() {
        std::lock_guard<std::mutex> lock(mtx_);
        finished_ = true;
        cv_.notify_all();
    }

private:
    std::queue<T> queue_;
    std::mutex mtx_;
    std::condition_variable cv_;
    bool finished_ = false;
};



// ---------- 并行处理主循环 ----------
template<typename First, typename Func>
void parallel_process(std::istream& in, size_t num_threads, First first, Func func) {
    TaskQueue<std::string> tasks;

    // 启动线程池
    std::vector<std::thread> workers;
    for (size_t i = 0; i < num_threads; ++i) {
        workers.emplace_back([&tasks, func]() {
            std::string chunk;
            while (tasks.pop(chunk)) {
                func(chunk); 
            }
        });
    }

    // 读取循环：保证每个块以换行结尾
    const size_t BLOCK_SIZE = 1 << 20;  // 1 MB
    std::vector<char> buffer(BLOCK_SIZE);
    std::string leftover;                // 上一个块末尾不完整的行
    bool firstDone = false;

    while (in.read(buffer.data(), buffer.size()) || in.gcount() > 0) {
        size_t bytes = in.gcount();
        std::string block(leftover + std::string(buffer.data(), bytes));

        // 定位最后一个换行符
        size_t last_nl = block.find_last_of('\n');
        if (last_nl == std::string::npos) {
            // 整个块没有换行，全部留到下一轮
            leftover = block;
            continue;
        }

        // 将完整部分（包括最后的换行）作为任务提交
        std::string complete = block.substr(0, last_nl + 1);
        leftover = block.substr(last_nl + 1);   // 剩余不完整行留到下次

        if (!complete.empty()) {
            if (!firstDone) {
                first(complete);
                firstDone = true;
            }
            tasks.push(std::move(complete));
        }
    }

    // 处理最后遗留的不完整行（可能没有换行符）
    if (!leftover.empty()) {
        // 如果最后一部分没有换行，也当作一个块处理
        // 因为数据流已结束，不会再有后续行拼接
        if (!firstDone) {
            first(leftover);
            firstDone = true;
        }
        tasks.push(std::move(leftover));
    }

    tasks.finish();
    for (auto& t : workers) {
        t.join();
    }
}



// ---------- 并行处理主循环 ----------
template<typename First, typename Func>
void parallel_process(const std::string& fname, size_t num_threads, First first, Func func) {
    std::istream* in = nullptr;
    std::unique_ptr<std::ifstream> file_ptr;

    if (fname == "-") {
        in = &std::cin;
    } else {
        file_ptr = std::make_unique<std::ifstream>(fname);
        if (!file_ptr->is_open()) {
            throw std::runtime_error("Cannot open file: " + fname);
        }
        in = file_ptr.get();
    }
    parallel_process(*in, num_threads, first, func);
}