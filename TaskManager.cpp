#include "TaskManager.h"

using json = nlohmann::json;

TaskManager* TaskManager::instance{nullptr};
std::mutex TaskManager::mutex;

TaskManager::~TaskManager() {
    shutdown();
}

TaskManager* TaskManager::GetInstance()
{
    std::lock_guard<std::mutex> lock(mutex);
    if (instance == nullptr)
    {
        instance = new TaskManager();
    }
    return instance;
}

void TaskManager::Start(){
    std::cout << "Task Manager starting!" << std::endl;

    if(mTaskThread.joinable()){
        std::cout << "Task Manager already running" << std::endl;
        return;
    }

    mTaskQueue.clear();
    //mTaskThread = std::thread(&TaskManager::run, this);
}

void TaskManager::QueueTask(Task task){
    // if(!mRunning){
    //     std::cout << "Task Manage is shut down. Please call start again" << std::endl;
    //     return;
    // }

    std::cout << "Queuing task " << task.type << std::endl;
    mTaskQueue.push_back(task);
}

void TaskManager::run(){
    std::cout << "Thread starting!" << std::endl;

    mRunning = true;
    while(mRunning){
        if(mTaskQueue.size() == 0) continue;

        execute(mTaskQueue.front());
        mTaskQueue.pop_front();
    }
    std::cout << "Thread shutting down!" << std::endl;
}

void TaskManager::execute(Task task){
    switch (task.type){
        case EXIT_THREAD:
            shutdown();
            break;
        case SEND_NODES:
            sendNodes(task.data, task.msg);
            break;
        case SEND_ROUTE:
            sendRoute(task.data, task.msg);
    }
}

void TaskManager::shutdown(){
    if(mTaskThread.joinable()){
        mRunning = false;
        mTaskThread.detach();
    }
}

void TaskManager::sendNodes(const std::vector<TA*>& nodes, const std::string msg){
    json j_array;
	json j;

	for(auto &item : nodes){
		json t;
		t["id"] = item->id;
		t["x"] = item->point.pos.lat;
		t["y"] = item->point.pos.lon;
		j_array.push_back(t);
		
	}

	j["type"] = "ta";
	j["list"] = j_array;

    std::cout << "Will send " << j.dump() << std::endl;

	const std::string str = j.dump();
	char * cstr = new char [str.length()+1];
  	std::strcpy (cstr, str.c_str());

    //send_message(cstr);
}

void TaskManager::sendRoute(const std::vector<TA*>& walk, const std::string msg){

    json j_array;
	json j;

    for(auto& ta : walk){
        j_array.push_back(ta->id);
    }

    j["dataType"] = "route";
	j["type"] = "ta";
	j["list"] = j_array;

    std::cout << "Will send " << j.dump() << std::endl;

	const std::string str = j.dump();
	char * cstr = new char [str.length()+1];
  	std::strcpy (cstr, str.c_str());

    //send_message(cstr);
}
