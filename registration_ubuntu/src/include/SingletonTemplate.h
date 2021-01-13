#pragma once

#include <memory>
#include <mutex>

namespace CloudReg
{

template<class T>
class rdSingleton
{
public:
    static T *instance();

private:
    rdSingleton()
    {
    }
    ~rdSingleton()
    {
    }

    static std::mutex lock_;
    static std::unique_ptr<T> instance_ptr_;
};

template<class T>
std::mutex rdSingleton<T>::lock_;

template<class T>
std::unique_ptr<T> rdSingleton<T>::instance_ptr_ = nullptr;

template<class T>
T *rdSingleton<T>::instance()
{
    if (nullptr == instance_ptr_)
    {
        std::lock_guard<std::mutex> guard(lock_);

        if (nullptr == instance_ptr_)
        {
            instance_ptr_.reset(new T());
        }
    }

    return instance_ptr_.get();
}

template <class T>
class WrapperSingleton
{
public:
    static T *instance(T *instance = nullptr)
    {
        static T *s_instance_ = nullptr;

        if (instance)
        {
            s_instance_ = instance;
        }

        return s_instance_;
    }
private:
    WrapperSingleton();
    ~WrapperSingleton ();
    WrapperSingleton(const WrapperSingleton &); // do not allow copy ctor.
    WrapperSingleton &operator=(const WrapperSingleton
                                &); // do not allow copy assign.
};

} 




