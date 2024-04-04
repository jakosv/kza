#include <chrono>

class Timer {
public:
    inline void start()
    {
        begin = std::chrono::steady_clock::now();
    }

    inline void stop()
    {
        end = std::chrono::steady_clock::now();
    }

    inline long elapsed_time()
    {
        auto diff = end - begin;
        return std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    }

    inline void print_elapsed(const char *msg)
    {
        long elapsed = this->elapsed_time();
        std::cout << msg << elapsed << "ms" << std::endl;
    }

private:
    std::chrono::steady_clock::time_point begin, end;
};
