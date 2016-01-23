/* 定义用户数据(userdata)结构体 */
struct userdata {
  double R; // 核半径
  double b; // 碰撞参量
  double x; // 场点横坐标
  double y; // 场点纵坐标
  double tau; // 固有时tau
  double Y0;  // 初始快度Y0
  double a;   // 参数a
  char flag;  // 核标记，'+'表示沿z轴正方向的核(即左核)，'-'表示沿z轴负方向的核(即右核)
  char type;  // 被积区域类型，'p'表示参与者，'s'表示旁观者
  double min[3]; // 积分下限，分别为x, y, Y的下限
  double max[3]; // 积分上限，分别为x, y, Y的下限
};
