#include <thread>
#include <deque>
#include <iostream>
#include <mutex>
#include "List.h"
#include "json.hpp"

#ifndef TASKMANAGER_H
#define TASKMANAGER_H

#define EXIT_THREAD -1
#define SEND_NODES 1
#define SEND_ROUTE 2

struct Task{
    int type;
    std::vector<TA*> data;
    std::string msg;
};

class TaskManager{
    public:
        //TaskManager should not be cloneable
        TaskManager(TaskManager &other) = delete;
        //TaskManager should not be assignable
        void operator=(const TaskManager &) = delete;

        static TaskManager *GetInstance();
        void Start();
        void QueueTask(Task);
        
    protected:
        TaskManager() = default;
        ~TaskManager();
    private:
        void run();
        void execute(Task);
        void shutdown();

        void sendNodes(const std::vector<TA*>&, const std::string);
        void sendRoute(const std::vector<TA*>&, const std::string);
    private:
        static TaskManager *instance;
        static std::mutex mutex;
        std::thread mTaskThread;
        std::deque<Task> mTaskQueue;
        bool mRunning = false;

        

};

#endif