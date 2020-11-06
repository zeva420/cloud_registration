OS  Win10-64
IDE VS2017 直接安装
三方库：
cmake 3.19.0 直接安装
PCL-1.8.1-AllInOne-msvc2017-win64.exe 直接安装
系统环境变量Path中增加如下配置
%PCL_ROOT%\bin
%PCL_ROOT%\3rdParty\VTK\bin
%PCL_ROOT%\3rdParty\FLANN\bin
%OPENNI2_REDIST64%
%PCL_ROOT%\3rdParty\Qhull\bin
%PCL_ROOT%\3rdParty\OpenNI2\Tools

g2o 请用提供的源代码包屏蔽了部分编译选项，eigen请用pcl安装时候携带的版本必须保持一致
1 修改eigen路径 CMakeList.txt中的 set(EIGEN3_INCLUDE_DIR "C:/Program Files/PCL 1.8.1/3rdParty/Eigen/eigen3") 
2 命令行进入g2o 解压目录 新建build目录后，cmake  -Thost=x64 ../ 生成VS工程
3 打开ALL_BUILD.vcxproj 生成release版本库 平台新建选择x64 如果链接器-》命令行-》其他选项”里面有/machine:X86直接删除



