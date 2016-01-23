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

/* 积分参数结构体 */
struct intargu {
  int ndim; // 积分维数
  int ncomp; // 被积函数分量数
  int nvec;  // 每次调用被积函数的最大计算点数， 一般设置为1
  double epsrel; // 相对误差限
  double epsabs; // 绝对误差限
  int flags; // 积分控制变量，一般用来设置输出信息的多少, 示例中的值为2
  int seed; // 随机种子， 示例中的值为0
  int mineval; // 被积函数被计算的最少次数，示例中的值为0
  int maxeval; // 被积函数被计算的(近似)最大次数，示例中的值为50000
  
  int nstart;  // [V]the number of integrand evaluations per iteration to start with. demo's value is 1000
  int nincrease; // [V]the increase in the number of integrand evaluations per iteration. demo's value is 500
  int nbatch; // [V]the batch size for sampling. demo's value is 1000
  int gridno; // [V]the slot in the internal grid table. demo's value is 0
  
  int nnew; // [S]the number of new integrand evaluations in each subdivision. demo's value is 1000
  int nmin; // [S]the minimum number of samples a former pass must contribute
            // to a subregion to be considered in that region’s compound integral value. demo's value is 2.
  double flatness; // [S]the parameter p in Eq. (1). demo's value is 25.

  int key1; // [D]determines sampling in the partitioning phase. demo's value is 47
  int key2; // [D]determines sampling in the final integration phase. demo's value is 1
  int key3; // [D]sets the strategy for the refinement phase. demo's value is 1
  int maxpass; // [D]controls the thoroughness of the partitioning phase. demo's value is 5.
  double border; // [D]the width of the border of the integration region. demo's value is 0.
  double maxchisq; // [D]the maximum χ 2 value a single subregion is allowed to have in the final integration phase. demo's value is 10.
  double mindeviation; // [D]a bound, given as the fraction of the requested error of the entire integral. demo's value is 0.25.
  int ngiven; // [D]the number of points in the xgiven array. demo's 0
  int ldxgiven; // [D]the leading dimension of xgiven, i.e. the offset between one point and the next in memory. demo's value is NDIM
  int nextra; // [D]the maximum number of extra points the peak-finder subroutine will return. demo's value is 0.

  int key; // [C]chooses the basic integration rule. demo's value is 0
  
  char *statefile; // 指定的状态文件名，示例中的值为 NULL
  void *spin; // the ‘spinning cores’ pointer, demo's value is NULL
  
};
