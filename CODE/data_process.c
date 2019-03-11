/** ###################################################################
**     Filename  : data_process.C
**     Project   : BalanceCar
**     Processor : MCF52255CAF80
**     Compiler  : CodeWarrior MCF C Compiler
**     Date/Time : 2013/9/26, 14:51
**     Contents  :
**         User source code
**
** ###################################################################*/

/* MODULE data_process */

#include "Cpu.h"
#include "Events.h"
#include "AS1.h"
#include "PWM_left.h"
#include "PWM_right.h"
#include "TI1MS.h"
#include "AD_ZOUT.h"
/* Including shared modules, which are used for whole project */
#include "PE_Types.h"
#include "PE_Error.h"
#include "PE_Const.h"
#include "IO_Map.h"
#include "data_pro.h"
////////////////////////////////////////////////////////////////////////

void GetInputVoltageAverage(void);              //函数原型声明
void Kalman_Filter(float angle_m,float gyro_m);
void AngleCalculate(void);
void SetMotorVoltage(float fLeftVoltage, float fRightVoltage);
void MotorSpeedOut(void);
Seria_data_process(float params,int params_id);


//////////////////////////////////////////////////////////////////////////

/***************************************/
/*查表数组 用于查-1~1的反三角函数asin()*/
/***************************************/
const float  Asin_to_Angle[201]=
{ 
-90.000000,-81.890386,-78.521659,-75.930132,-73.739795,-71.805128,-70.051556,-68.434815,-66.926082,-65.505352,
-64.158067,-62.873247,-61.642363,-60.458639,-59.316583,-58.211669,-57.140120,-56.098738,-55.084794,-54.095931,
-53.130102,-52.185511,-51.260575,-50.353889,-49.464198,-48.590378,-47.731416,-46.886394,-46.054480,-45.234915,
-44.427004,-43.630109,-42.843643,-42.067065,-41.299873,-40.541602,-39.791819,-39.050123,-38.316134,-37.589503,
-36.869898,-36.157008,-35.450543,-34.750226,-34.055798,-33.367013,-32.683639,-32.005455,-31.332251,-30.663830,
-30.000000,-29.340582,-28.685402,-28.034297,-27.387108,-26.743684,-26.103881,-25.467560,-24.834587,-24.204835,
-23.578178,-22.954499,-22.333683,-21.715617,-21.100196,-20.487315,-19.876874,-19.268775,-18.662925,-18.059230,
-17.457603,-16.857956,-16.260205,-15.664267,-15.070062,-14.477512,-13.886540,-13.297072,-12.709033,-12.122352,
-11.536959,-10.952784,-10.369760,-9.7878190,-9.2068960,-8.6269270,-8.0478460,-7.4695920,-6.8921030,-6.3153160,
-5.7391700,-5.1636070,-4.5885660,-4.0139870,-3.4398130,-2.8659840,-2.2924430,-1.7191310,-1.1459920,-0.5729670,
0.00000000,
0.57296700,1.14599200,1.71913100,2.29244300,2.86598400,3.43981300,4.01398700,4.58856600,5.16360700,5.73917000,
6.31531600,6.89210300,7.46959200,8.04784600,8.62692700,9.20689600,9.78781900,10.3697600,10.9527840,11.5369590,
12.1223520,12.7090330,13.2970720,13.8865400,14.4775120,15.0700620,15.6642670,16.2602050,16.8579560,17.4576030,
18.0592300,18.6629250,19.2687750,19.8768740,20.4873150,21.1001960,21.7156170,22.3336830,22.9544990,23.5781780,
24.2048350,24.8345870,25.4675600,26.1038810,26.7436840,27.3871080,28.0342970,28.6854020,29.3405820,30.0000000,
30.6638300,31.3322510,32.0054550,32.6836390,33.3670130,34.0557980,34.7502260,35.4505430,36.1570080,36.8698980,
37.5895030,38.3161340,39.0501230,39.7918190,40.5416020,41.2998730,42.0670650,42.8436430,43.6301090,44.4270040,
45.2349150,46.0544800,46.8863940,47.7314160,48.5903780,49.4641980,50.3538890,51.2605750,52.1855110,53.1301020,
54.0959310,55.0847940,56.0987380,57.1401200,58.2116690,59.3165830,60.4586390,61.6423630,62.8732470,64.1580670,
65.5053520,66.9260820,68.4348150,70.0515560,71.8051280,73.7397950,75.9301320,78.5216590,81.8903860,90.0000000
};

//////////////////////////////////////////////////////////////////////////

unsigned short g_nADGravity[10]={0};          //存放采集到的直接AD值
unsigned short g_nADGyroscope[10]={0};

unsigned int g_nGravityVoltage=0;         //存储AD的平均值
unsigned int g_nGyroVoltage=0;

volatile float gravity_offset= 2053.0;      //z轴零偏值
volatile float gryoscope_offset=1698.0;     //陀螺仪零偏值
volatile float gravity_angle_ratio=-0.08923;  //单位换算比例系数
volatile float gryoscope_angle_ratio=0.24;
volatile float left_dead_val=0;           //死区电压
volatile float right_dead_val=0;
volatile float angle_control_p=0;
volatile float angle_control_d=0;
volatile float g_fGyroscopeAngleSpeed=0;    //转化后的车体实际物理加速度
volatile float g_fCarAngle=0;               //转化后的实际车体物理倾斜角

volatile float Angle, Angle_dot; 		//外部需要引用的变量

word speed=32768;
////////////////////////////////////////////

volatile float Q_angle=0.001, Q_gyro=0.003, R_angle=0.67, dt=0.005;   	//dt的取值为kalman滤波器采样时间;		
volatile float P[2][2] ={{ 1, 0 },{ 0, 1 }};  ///////////////////
volatile char  C_0 = 1;
volatile float E;  ///////////////////////////
volatile float q_bias;/////////////////////
volatile float Angle_err;	//////////////////////
volatile float PCt_0, PCt_1;  //////////////////
volatile float K_0, K_1;//////////////////////
volatile float t_0, t_1;///////////////////////
volatile float Pdot[4] ={0,0,0,0};	////////////////
////////////////////////////////////////////////////////
bool g_nRunning_flag=1;
volatile float g_fAngleControlOut=0;        //平衡控制理想电压值

volatile float g_fLeftMotorOut=0;           //上面三个理想电压值的叠加量
volatile float g_fRightMotorOut=0;          //上面三个理想电压值的叠加量



void SampleInputVoltage()               //直接AD采样
{
	byte i=0;
	for(i=0;i<10;i++)
	{
		AD_ZOUT_Measure(1);
		AD_ZOUT_GetChanValue16(0,&(g_nADGravity[i]));
	  	AD_ZOUT_GetChanValue16(1,&(g_nADGyroscope[i])); 
	}
	

}


void GetInputVoltageAverage(void)                   //求AD采样的平均值 
{
	byte i=0;
	unsigned long AD_Gravity_sum=0,AD_Gyroscope_sum=0;
	for(i=0;i<10;i++)
	{
		AD_Gravity_sum+=(unsigned long)g_nADGravity[i];
		AD_Gyroscope_sum+=(unsigned long)g_nADGyroscope[i];
	}
	g_nGravityVoltage=(AD_Gravity_sum>>4)/10;
	g_nGyroVoltage=(AD_Gyroscope_sum>>4)/10;
}


void Kalman_Filter(float angle_m,float gyro_m)			//gyro_m:gyro_measure
{
	Angle+=(gyro_m-q_bias) * dt;    //先验估计
	
	Pdot[0]=Q_angle - P[0][1] - P[1][0];// Pk-' 先验估计误差协方差的微分
	Pdot[1]=- P[1][1];
	Pdot[2]=- P[1][1];
	Pdot[3]=Q_gyro;
	
	P[0][0] += Pdot[0] * dt;        // Pk- 先验估计误差协方差微分的积分 = 先验估计误差协方差
	P[0][1] += Pdot[1] * dt;
	P[1][0] += Pdot[2] * dt;
	P[1][1] += Pdot[3] * dt;
		
	Angle_err = angle_m - Angle;    //zk-先验估计
		
	PCt_0 = C_0 * P[0][0];
	PCt_1 = C_0 * P[1][0];
	
	E = R_angle + C_0 * PCt_0;
	
	K_0 = PCt_0 / E;                //Kk
	K_1 = PCt_1 / E;
	
	t_0 = PCt_0;
	t_1 = C_0 * P[0][1];

	P[0][0] -= K_0 * t_0;           //后验估计误差协方差
	P[0][1] -= K_0 * t_1;
	P[1][0] -= K_1 * t_0;
	P[1][1] -= K_1 * t_1;
		
	Angle	+= K_0 * Angle_err;       //后验估计	
	q_bias+= K_1 * Angle_err;       //后验估计
	Angle_dot = gyro_m-q_bias;      //输出值（后验估计）的微分 = 角速度
}


#define GRAVITY_OFFSET           gravity_offset       //z轴零偏值
#define GYROSCOPE_OFFSET         gryoscope_offset     //陀螺仪零偏值
#define GRAVITY_ANGLE_RATIO      gravity_angle_ratio  //单位换算比例系数
#define GYROSCOPE_ANGLE_RATIO    gryoscope_angle_ratio
#define ANGLE_INT                29.6//-40.2698    
 
void AngleCalculate(void)                          //角度计算函数
{
	float Gyro_m = 0;
    float Angle_m = 0;
    float Ang_Value;
    
	//==========加速度器处理部分==========//
	int Ang_i;
	Ang_Value=(float)((int)g_nGravityVoltage-(int)gravity_offset);
	Ang_Value=Ang_Value * (3300/4096);           //获得电压值
	Ang_Value = Ang_Value / 7.60;                //求出 Z 轴加速度 求出sin值 然后放大100倍 PS.传感器800mv/g 实测760mv/g
	Ang_i=(int)Ang_Value;
	if(Ang_i>100) Ang_i=100;                     //限幅
    if(Ang_i<-100) Ang_i=-100;
    Angle_m = Asin_to_Angle[(uint8_t)(Ang_i + 100)];       //计算角度
    
    //==========陀螺仪处理部分==========//   
    Ang_Value = (float)((int)g_nGyroVoltage - (int)gryoscope_offset); //减去零偏置值
    Ang_Value = Ang_Value * (3300/4096);                              //获得电压值
    Gyro_m = Ang_Value / (670*5.1);                   //计算角角速度,陀螺仪灵敏度：0.67mv/deg/sec

    //======= 卡尔曼滤波器 ========//
    Kalman_Filter(Angle_m, Gyro_m);
        
	g_fGyroscopeAngleSpeed=Angle_dot;
	g_fCarAngle=Angle;
}


#define ANGLE_CONTROL_P  angle_control_p
#define ANGLE_CONTROL_D  angle_control_d

void AngleControl(void)        //角度控制函数 使用的是PD控制  
{
	float fValue;
	GetInputVoltageAverage();
	AngleCalculate();
    fValue =(-(g_fCarAngle-ANGLE_INT))* ANGLE_CONTROL_P +(-g_fGyroscopeAngleSpeed) * ANGLE_CONTROL_D;
	g_fAngleControlOut = fValue;	
}


void SetMotorVoltage(float fLeftVoltage, float fRightVoltage)
{
	float Lduty=32768;
	float Rduty=32768;
	Lduty=32768+fLeftVoltage*4550;
	Rduty=32768+fRightVoltage*4550;
	if(g_nRunning_flag==1)
	{
	if(Lduty>65535)
		Lduty=65535;
	if(Lduty<0)
		Lduty=0;
	if(Rduty>65535)
		Rduty=65535;
	if(Rduty<0)
		Rduty=0;
		PWM_left_SetRatio16((word)Lduty);
		PWM_right_SetRatio16((word)Rduty);
	}
	else
	{
		PWM_left_SetRatio16(32767);
		PWM_right_SetRatio16(32767);
	}
}

#define LEFT_OUT_DEAD_VAL  left_dead_val
#define RIGHT_OUT_DEAD_VAL  right_dead_val
void MotorSpeedOut(void) 
{
	float fLeftVal, fRightVal;
	
	fLeftVal = g_fLeftMotorOut;
	fRightVal = g_fRightMotorOut;
	if(fLeftVal > 0) 			fLeftVal += LEFT_OUT_DEAD_VAL;
	else if(fLeftVal < 0) 		fLeftVal -= LEFT_OUT_DEAD_VAL;
	
	if(fRightVal > 0)			fRightVal += RIGHT_OUT_DEAD_VAL;
	else if(fRightVal < 0)		fRightVal -= RIGHT_OUT_DEAD_VAL;			
	SetMotorVoltage(fLeftVal, fRightVal);
}

#define MOTOR_OUT_MAX  7.2
#define MOTOR_OUT_MIN  -7.2

void MotorOutput(void) 
{
	float fLeft, fRight;

	fLeft = g_fAngleControlOut;
	fRight = g_fAngleControlOut;
	
	if(fLeft > MOTOR_OUT_MAX)		fLeft = MOTOR_OUT_MAX;
	if(fLeft < MOTOR_OUT_MIN)		fLeft = MOTOR_OUT_MIN;
	if(fRight > MOTOR_OUT_MAX)		fRight = MOTOR_OUT_MAX;
	if(fRight < MOTOR_OUT_MIN)		fRight = MOTOR_OUT_MIN;
		
	g_fLeftMotorOut = fLeft;
	g_fRightMotorOut = fRight;
	MotorSpeedOut();
}








 Seria_data_process(float params,int params_id)  //串口接收
{
	/////按顺序来/////////
	if(params_id == 0)
		angle_control_p = params;
	if(params_id == 1)
		angle_control_d = params;
}


















/* END data_process */
