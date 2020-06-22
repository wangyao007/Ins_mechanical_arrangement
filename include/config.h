#ifndef CONFIG_H
#define CONFIG_H

#include"opencv2/core/core.hpp"
#include<string>
#include <memory>
#include<iostream>
class Config
{
  private:
  static std::shared_ptr<Config> config_;
  cv::FileStorage file_;
  Config(){}
  public:
  ~Config();
  static void setParameterFile(const std::string& filename);
  template< typename T >   //模板函数，用于传递各类初始值
  static T get( const std::string& key )
    {
        return T( Config::config_->file_[key] );
    }
};  
#endif
