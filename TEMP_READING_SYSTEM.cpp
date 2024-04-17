#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstdlib>
using namespace std;
int MODULE_CALLS();
int MODULE1();
int MODULE2();
int MODULE3();
int MODULE4();
int MODULE5();
int MODULE6();
int MODULE7();
int temp_conversions();
int cel_fahren();
int cel_kelvin();
int fahren_kelvin();
int solid_formation_table(double,double,double,double,double);
int Temp_of_two_liquids();
int Phases_of_water();
int Ice_and_Water();
int Water_and_Steam();
int Temperature_of_cooling_body();
int Temperature_chart(double,double,double,double);
int Estimation_of_the_time();
int converting_into_newscale();
class metal_bar;    
/*class TEMPERATURE_READING_SYSTEM();
{
 public:
 TEMPERATURE_READING_SYSTEM(){}*/
int main()
{
    MODULE_CALLS();
}
  int MODULE_CALLS()
   {
    int s;
    do
    {
     std::system("clear");
     cout<<"\n\n\t\t\t\t\tTEMPERATURE READING SYSTEM\n\n\t\t\t";
     int choose_option;
     cout<<"1.NEWTON'S LAW OF COOLING APPLICATIONS\n\n\t\t\t";
     cout<<"2.JUNCTION TEMPERATURE BETWEEN TWO METAL RODS\n\n\t\t\t";
     cout<<"3.MIXTURE TEMPERATURE OF MISCIBLE LIQUIDS AND SOLIDS\n\n\t\t\t";
     cout<<"4.TEMPERATURE MAINTAINED TO FREEZE A SOLID UNDER ITS FUSION TEMPERATURE WITHIN CERTAIN TIME\n\n\t\t\t";
     cout<<"5.ESTIMATION THE TEMPERATURE THAT IS TO BE MAINTAINDE DURING THERMAL EXPANSION OF A METAL ROD\n\n\t\t\t";
     cout<<"6.TEMPERATURE FORECAST OF AREGION BASED ON PRESSURE VARIATIONS AT SEA LEVEL\n\n\t\t\t";
     cout<<"7.TEMPERATURE INTERCONVERSIONS\n\n\t\t\t";
     cout<<"\n\n\t\tEnter the option from (1-7):";
     cin>>choose_option;
     switch(choose_option)
     {
      case 1:
      {
        MODULE1();
       break;
      }
      case 2:
      {
        MODULE2();
       break;
      }
      case 3:
      {
        MODULE3();
       break;
      }
      case 4:
      {
        MODULE4();
       break;
      }
      case 5:
      {
        MODULE5();
       break;
      }
      case 6:
      {
        MODULE6();
       break;
      }
      case 7:
      {
        MODULE7();
       break;
      }
      default:
      {
       cout<<"\n\t\t\tDo Enter the right choice";
       break;
      }
     }  
     cout<<"\n\n\t\t\tTo go to Home page Press 1:";
     cin>>s;
     }while(s==1);
  }
int MODULE1()
{
 int x=0;
    do
    {
    std::system("clear");
    int choice=0;
    cout<<endl<<"\t\t\tNEWTON'S LAW OF COOLING APPLICATIONS"<<endl<<endl;
    cout<<"\t\t1.To Find Temperature of body at particular time while cooling"<<endl;
    cout<<"\t\t2.Estimation of the time of death of a body based on Temperature readings[Application]"<<endl;
    cout<<endl<<"\t\tEnter the correct choice:";
    cin>>choice;
    cout<<endl;
    switch(choice)
    {
        case 1:
        {
         cout<<"\t\tTo Find the Temperature of a Body at particular time while cooling"<<endl; 
         Temperature_of_cooling_body();
         break;
        }
        case 2:
        {
            cout<<"\t\tEstimation of Time of Death"<<endl;
            Estimation_of_the_time();
            break;
        }
        default:
        {
         cout<<"\t\tChoose the correct choice i.e from (1-2)"<<endl;
         break;
        }
    }
    cout<<"\n\t\tTo return to Main Page press 1:";
    cin>>x;
    }while(x==1);
}
int Temperature_of_cooling_body()
{
    int y=0;
    do
    {
     std::system("clear");
    //This function gives the temperature of a cooling body at 
    //a particular time[considering the time from when it is left to cool down]
    double T1,T2,T;
    //T1 is the initial temperature
    //T2 is the temperature at mentioned time(observation)
    //T is the surrounding temperature
    //For Higher temperatures we use Stefan's Law
    //Newton's laws of coolin is applicable in day to day activities
    //IMPORTANT NOTE:The TIME provided should be from the time when it is left to cool down[t=0]
    double Time1,Time2;
    //Time1 is the observation time
    //Time2 is the given time where temperature is to be found out
    cout<<"\n\t\t\tTo Find the Temperature of a Body at particular time while cooling"<<endl<<endl; 
    cout<<"\t\tProvide the following details from observations:"<<endl<<"\n\t\t";
    cout<<"Initial Temperature of the body('C)=";
    cin>>T1;
    cout<<"\n\t\tSurroundings Temperature('C)=";
    cin>>T;
    if(T1<=T)
    {
    cout<<"\n\t\tSince Surrounding Temperature is greater/Equal than the body Temperature the body doesn't\n\t\t";
    cout<<"cool and Temperature of body remains same i.e "<<T1;
    }
    else
    {    
    cout<<"\n\t\tInitial obsevartions:"<<endl<<"\n\t\t";
    cout<<"Time lapsed(during observation)(min)=";
    cin>>Time1;
    cout<<"\n\t\tTemperature('C) after Observation Time=";
    cin>>T2;
    double T3;
    cout<<endl<<"\t\tTime at which temperature is required(min)=";
    cin>>Time2;
    double k,c;
    //k and c constants
    k=(2*(T1-T2))/(Time1*(T1+T2-(2*T)));
    c=(k*Time2)/2;
    T3=((T2*(1-c))+(2*c*T))/(1+c);
    cout<<"\n\t\tThe Temperature('C) of Body  after "<<Time2<<"min is "<<T3<<endl;
    int x=0;
    cout<<"\n\t\tTo print the record chart of Temperatures at succesive intervals of time press 1:";
    cin>>x;
    cout<<endl<<"\t\t";
    if(x==1)
    {   
        double interval;
        cout<<"\t\tTime interval(min)=";
        cin>>interval;
        Temperature_chart(T1,k,T,interval);
    }
    }
    cout<<endl<<"\n\t\tTo repeat the process for different values press 1:";
    cin>>y;
    }while(y==1); 
}
int Temperature_chart(double Temp1,double k,double T,double time)
{
   int x=0;
   do
   {
   double T2,i=0.0;
   cout<<"\t\tEnter till how many intervals the chart should display"<<endl<<"\t\t";
   int n;
   cin>>n;
   cout<<endl<<endl<<"\t\tTIME(min)\t\tTEMPERATURE('C)"<<endl;
   for(int j=0;j<n;j++)
   {
       T2=(Temp1*(2-(k*time))+(2*k*time*T))/(2+(k*time));
       i=i+time;
       Temp1=T2;
       cout<<"\t\t"<<i<<"\t\t\t"<<T2<<endl;
   }
   cout<<endl<<"\t\tTo produce more values in chart press 1 and increase the interval:";
   cin>>x;
   }while(x==1);
}
int Estimation_of_the_time()
{    
     int x;
     do
     {
       std::system("clear");
       cout<<"\n\t\t\tEstimation of Time of Death"<<endl;
    //This function is used mostly in forensic labs to estimate the time of death of a person[body]
       cout<<endl<<"\t\tProvide the following details"<<endl<<"\t\t";
       double temp1,temp2,temp3,T;
       //temp1 is the normal body temperature considering the person was clinically normal i.e 37('C)
       //temp2 is the first obsevation at time t1 after death
       //temp3 is the second obsevation at time t2 after death
       //T is the surrounding temperature
       double interval,time;
       //interval is t2-t1
       cout<<"\n\t\tInitial observation of temperature('C)=";
       cin>>temp2;
       cout<<"\t\tFinal observation of Temperature('C)[after time interval)=";
       cin>>temp3;
       cout<<"\t\tTime interval between two observations(min)=";
       cin>>interval;
       //time is the time of death before the arrival of forensic
       cout<<"\t\tSurrounding Temperature('C)=";
       cin>>T;
       double k;
       //k is constant
       if(temp2<=37&&temp2>temp3&&temp3<=37)
       {
       k=(2*(temp2-temp3))/(interval*(temp2+temp3-(2*T)));
       //37 is the normal body temperature
       time=(2*(37-temp2))/(k*(37+temp2-(2*T)));
       cout<<"\n\n\t\tEstimated time of the death is "<<time<<" min before your arrival"<<endl<<endl;
       }
       else
       {
       cout<<"\n\t\tThe values entered are not correct.Enter correct readings.";
       cout<<"\n\t\tNOTE:The body Temperature decreases after the person expires\n\n";
       }
       cout<<"\t\tTo continue for further computations press 1:";
       cin>>x;
     }while(x==1);
}
class metal_bar
{
public:
 double A;
 double L;
 double K;
 metal_bar(){}
 void insert(double area,double length,double conductivity)
 {
    A=area;
    L=length;
    K=conductivity;
 }
friend double junction_temp(metal_bar r1,metal_bar r2,double T1,double T2);
friend int MODULE2();
}; 
int MODULE2()
{
  int x=0;
    do
    {
    std::system("clear");
    cout<<endl;
    double Temp1=0;
    double Temp2=0;
    double area,length,conductivity;
    metal_bar rod[2];
    //Temp1 is the temperature at one end of rod1 and Temp2 is the temperature at aother end of rod2
    //The purpose is to calculate the junction temperture between two rods when joined 
    //considering the same rate
    cout<<"\t"<<"\t"<<"JUCTION TEMPERATURE BETWEEN TWO RODS"<<endl<<endl;
    cout<<"\t"<<"Provide the following details:"<<endl;
    cout<<"\t [ Some Thermal onductivity values-\n\t   Copper=386\n\t   Steel=46\n\t ]"<<endl<<endl;
    for(int i=0;i<2;i++)
    {
        cout<<"\t"<<"\t"<<"METAL ROD "<<i+1<<endl;
        cout<<"\t"<<"conductivity(J/sm'c)=";
        cin>>conductivity;
        cout<<endl<<"\t"<<"area of crossection(cm^2)=";
        cin>>area;
        cout<<endl<<"\t"<<"length(cm)=";
        cin>>length;
        cout<<endl<<endl;
        rod[i].insert(area,length,conductivity);
    }  
    cout<<"\t\tMention the temperatures maintained at free ends of the rods after joining them:"<<endl;
    cout<<"\n\t\tTemp1('C)=";
    cin>>Temp1;
    cout<<endl<<"\t\tTemp2('C)=";
    cin>>Temp2;
    cout<<endl;
    cout<<"\t"<<"The Junction Temperature of two metal rods is('C): ";
    junction_temp(rod[0],rod[1],Temp1,Temp2);
    cout<<"\n\t\tIf you want to continue for further computations press 1:";
    cin>>x;
    }while(x==1);
}
double junction_temp(metal_bar r1,metal_bar r2,double T1,double T2)
{
    double x=0.0;
    double y;
    x=((T2*r2.K*r2.A*r1.L)+(T1*r1.K*r1.A*r2.L));
    y=((r1.K*r1.A*r2.L)+(r2.K*r2.A*r1.L));
    double z=x/y;
    cout<<z<<endl;
    return 0;
}
int MODULE3()
{     
  int x=0;
    do
    {
    std::system("clear");
    int choice;
    cout<<"\n\t\t\tMIXTUTE TEMPERATURE OF MISCIBLE LIQUIDS AND SOLIDS"<<endl<<"\t\t";
    cout<<endl<<"\n\t\t1.Mixture Temperaure of two liquids"<<endl<<"\t\t";
    cout<<"2.Mixture Temperature of different phases of Water[ice,water,steam]"<<endl<<"\n\t\t";
    cout<<"Enter the choice:";
    cin>>choice;
    switch(choice)
    {
        case 1:
        {
            cout<<endl<<"\t\t Mixture Temperaure of two miscible liquids "<<endl;
            Temp_of_two_liquids();
            break;
        }
        case 2:
        {
            cout<<endl<<"\t\tMixture Temperature of combinations of ice,water,steam"<<endl;
            Phases_of_water();
            break;
        }
        default:
        {
            cout<<endl<<"\t\tChoose the correct option (1-2)";
            break;
        }
    }
    cout<<"\n\t\tTo return to Main page press 1:";
    cin>>x;
    cout<<endl;
    }while(x==1);
}
int Temp_of_two_liquids()
{   
    int x=0;
    do
    {
    double M1,M2,C1,C2,T1,T2,T;
    //M1 and M2 are the masses of liquids taken
    //C1 and C2 are the specific heats of the liquids
    //T1 and T2 are the individual Temperatures of the liquids
    cout<<endl<<"\t\tProvide the  Following details[T1>T2]:"<<endl;
    cout<<"\t\tLiquid 1:"<<endl;
    cout<<"\t\tMass(gm)=";
    cin>>M1;
    cout<<"\t\tTemperature('C)=";
    cin>>T1;
    cout<<"\t\tSpecific heat(cal/gm.'C)=";
    cin>>C1;
    cout<<endl<<"\t\tLiquid 2:"<<endl;
    cout<<"\t\tMass(gm)=";
    cin>>M2;
    cout<<"\t\tTemperature('C)=";
    cin>>T2;
    cout<<"\t\tSpecific heat(cal/gm.'C)=";
    cin>>C2;
    double k;
    //k is constant
    k=(M2*C2)/(M1*C1);
    double t;
    //t is common[mixture] Temperaure of two liquids
    t=(T1+(k*T2))/(1+k);
    cout<<endl<<"\t\tThe Mixture Temperature('C) is "<<t<<endl;
    cout<<"\n\t\tTo repeat the computations wuth different values press 1:";
    cin>>x;
    }while(x==1);
}
int Phases_of_water()
{
 int option,y=0;
 do
 {
 cout<<endl<<"\t\t1.Ice and Water mixture"<<endl;
 cout<<"\t\t2.Water and Steam mixture"<<endl;
 cout<<endl<<"\t\tChoose the Option:";
 cin>>option;
 switch(option)
 {
     case 1:
    { 
      cout<<"\n\t\tIce and Water Mixture\n\n\t\t";
      Ice_and_Water();
      break;
    }
    case 2:
    {
      cout<<"\n\t\tWater and Steam mixture\n\n\t\t";
      Water_and_Steam();
      break;
    }
    default:
    {
      cout<<"\t\tChoose the right option from (1-2)\n";
      break;
    }
 }
 cout<<"\n\t\tTo repeat the computaions press 1:";
 cin>>y;
 }while(y==1);
}
int Ice_and_Water()
{
    double L=80,S1=0.5,S2=1;
    //L is the latent heat(cal/gm) of Fusion
    //S1 is the specifc heat(cal/gm.'C) of ice
    //S2 is the specificheat of water
    double M1,M2,T1,T2,T;
    cout<<"Provide the details:\n\n\t\t";
    cout<<"Ice:\n\t\t";
    cout<<"Mass(gm)=";
    cin>>M1;
    cout<<"\t\tTemperature('C)=";
    cin>>T1;
    cout<<"\n\t\tWater:\n\t\t";
    cout<<"Mass(gm)=";
    cin>>M2;
    cout<<"\t\tTemperature('C)=";
    cin>>T2;
    double H1,H2,H;
    //H1 is the heat required to convert ice into water
    //H2 is the heat reqired by water to reach 0'C
    H2=(M2*S2*T2);
    H1=(M1*S1*(-T1))+(M1*L);
    if(H1>=H2)
    {
        double m;
        H=H1-H2;
        m=H/L;
        //The whole ice is not melted hence the Mixture Temperature is 0
        cout<<"\n\n\t\tThe whole ice is not melted hence,the Temperature of the mixute is 0('C)"<<endl;
        cout<<"\t\tMass of ice in the mixture is "<<m<<"gm"<<endl;
        cout<<"\t\tMass of water in the mixture is "<<(M2+(M1-m))<<" gm"<<endl;
    }
    else
    {
      double m,k;
      // is the mass of water remained at T2
      H=H2-H1;
      m=(H/H2)*M2;
      k=m/M1;
      T=(k*T2)/(1+k);
      cout<<"\n\n\t\tThe whole ice is melted and the Temperature of the mixture is "<<T<<"('C)"<<endl;
      cout<<"\t\tMass of ice in mixture is 0 gm"<<endl;
      cout<<"\t\tMass of water in the mixture is "<<(M1+M2)<<" gm"<<endl;
    }
}
int Water_and_Steam()
{
   double L=539,S2=1;
   //L is the latent heat(cal/gm) of vaporization
   //S2 is the specific heat of water
   double M1,M2,T1=100,T2,T;
   cout<<"Provide the details:\n\n\t\t";
   cout<<"Steam:\n\t\t";
   cout<<"Mass(gm)=";
   cin>>M1;
   cout<<"\t\tTemperature('C)=100";
   cout<<"\n\t\tWater:\n\t\t";
   cout<<"Mass(gm)=";
   cin>>M2;
   cout<<"\t\tTemperature('C)=";
   cin>>T2;
   double H1,H2;
   //H1 is the heat required to convert steam into water
   //H2 is the heat reqired by water to reach 100
   H2=(M2*S2*(T1-T2));
   H1=(M1*L);
   if(H1>=H2)
   {
    double m;
    m=(H2/L);
    //The whole steamis not condensed to water, hence the Mixture Temperature is 100('c)
    cout<<"\n\n\t\tThe whole steam is not condensed to water hence,the Temperature of the mixute is 100('C)"<<endl;
    cout<<"\t\tMass of steam in the mixture is "<<(M1-m)<<"gm"<<endl;
    cout<<"\t\tMass of water in the mixture is "<<(M2+m)<<" gm"<<endl;
   }
   else
   {
     double k;
     k=(M2/M1);
     T=((k*T2)+T1+L)/(1+k);
      cout<<"\n\n\t\tThe whole steam is condensed and the Temperature of the mixture is "<<T<<"('C)"<<endl;
      cout<<"\t\tMass of steam in mixture is 0 gm"<<endl;
      cout<<"\t\tMass of water in the mixture is "<<(M1+M2)<<" gm"<<endl;
   }
}
                                                                                        
int MODULE4()
{
  int z=0;
    do
    {
     std::system("clear");
    double K;
    double p,L;
    //K is the Themal conductivity of solid
    //p is the density of the solid beig formed
    //L is the latent heat of fusion of solid
    cout<<"\n\t\tTEMPERATURE TO BE MAINTAINED TO FREEZE A SOLID UNDER ITS FUSION POINT WTHIN CERTAIN TIME"<<endl;
    cout<<"\n\t\tProvide the following details:"<<endl;
    cout<<"\t\t [For Ice-\n\t\t  Latent heat=80\n\t\t  Density=0.92\n\t\t  Thermal conductivity=0.004\n\t\t ]"<<endl;
    cout<<"\n\t\tThermal conductivity of solid(ergs/cm.s.C)=";
    cin>>K;
    cout<<endl<<"\t\tDensity of solid(gm.cm^-3)=";
    cin>>p;
    cout<<endl<<"\t\tLatent Heat of Fusion(cal/gm)=";
    cin>>L;
    double Temp;
    //Temp is the external Tempertaure to be maintained 
    //Make sure you maintain the other end at fusion Temperature
    double constant=(p*L)/(2*K);
    cout<<endl<<"\t\tEnter the Fusion Temp of the solid('C)=";
    double Ftemp,H1;
    cin>>Ftemp;
    cout<<endl<<"\t\tEnter the desired Height of Block and the time at which the solid is to be freezed"<<endl;
    cout<<"\t\tInitial height(cm)=";
    cin>>H1;
    cout<<endl<<"\t\tDesired Height(cm)=";
    double H2,T;
    cin>>H2;
    cout<<endl<<"\t\tTime at which the solid is to be freezed(hrs)=";
    cin>>T;
    T=T*3600;
    Temp=Ftemp-(constant*((H2*H2)-(H1*H1)))/T;
    cout<<endl<<"\t\tThe external temperature('C)to be maintained is:"<<Temp<<endl;
    cout<<endl<<"\t\tTo display the Height status Report press 1"<<endl<<"\t\t";
    int x=0;
    cin>>x;
    if(x==1)
    { 
        int y;
        do
        {
        solid_formation_table(constant,T,Temp,H2,H1);
        cout<<"\t\t To find at different intervals press 1"<<endl<<"\t\t";
        cin>>y;
        }while(y==1);
    }
   cout<<"\t\t If want to continue for further computations press 1"<<endl<<"\t\t";
   cin>>z;
    }while(z==1);
}
int solid_formation_table(double constant,double time,double Temp,double H2,double H1)
{
  double interval,t=0.0,p=0.0;
  cout<<endl<<"\t\tEnter the time interval(hrs):";
  cin>>interval;
  cout<<endl<<"\t\tTIME(hrs)\t\tHEIGHT(cm)";
  double height=H1;
  while(height<H2&&t<time)
  {
    p=(H1*H1)+(-(Temp*t))/constant;
    height=sqrt(p);
    cout<<endl<<"\t\t"<<(t/3600)<<"\t\t\t"<<height<<endl;;
    t=t/3600;
    t=t+interval;
    t=t*3600;
    height=height;
  }
  cout<<endl;
}
int MODULE5()
{                                                                                                        
  int x=0;
    do
    {
    std::system("clear");
    cout<<endl;
    cout<<"\t\t"<<"Estimating the Temperature that is maintained during Thermal expansion of a metal rod"<<endl;
    double L1,L2,T1,T2,alpha;
    //L1 is the intial length at time 0
    //L2 is the final length at time t
    //T1 is the surrounding temp
    //alpha is the thermal coefficient of the metal
    //T2 is the temperature that is to be maintained
    cout<<endl<<"\t"<<"Provide the following details:"<<endl;
    cout<<"\t[Some Thermal conductivity(*10^-6) values-"<<endl;
    cout<<"\t   Iron=12\n\t   Copper=17\n\t   Aluminium=24\n\t]\n\n";
    cout<<"\t"<<"Temperature of surroundings('C)=";
    cin>>T1;
    cout<<"\t"<<"Initial length(cm)=";
    cin>>L1;
    cout<<"\t"<<"Final lenght(cm)=";
    cin>>L2;
    cout<<"\t"<<"Thermal coefficient(*10^-6'C^-1)=";
    cin>>alpha;
    T2=T1+((L2-L1)/(L1*alpha))*pow(10,6);
    cout<<endl<<"\t"<<"The Temperature('C) to be maintained to get the desired final length is:";
    cout<<T2<<endl<<endl;
    cout<<"\t\tIf you want to do for further computations press 1:";
    cin>>x;
    }while(x==1);
}
//This function is mainly used to know the
//temperature reading due to pressure changes
//at sea level due to changes in weather
int MODULE6()
{
   int x;
 do
 {
 std::system("clear");
 cout<<endl;
 double h,d=1,g=9.8,T1;
 //h is the height of a region
 //d is the density of the region[ideal case is chosen in case of density]
 //g is the gravity[no variation with respect to height...h<<R/2],whre R is the radius of Earth
 cout<<"\t\t\tTEMPERATURE FORECAST OF A REGION BASED ON PRESSURE VARIATIONS AT SEA LEVEL"<<endl<<endl;
 cout<<"\t\tProvide the following details based on Observations:"<<endl<<endl;
 double Patm1=101325;
 //Patm1 is the atmospheric pressure at ideal conditions
 //1 atm=101325 pascal
 cout<<"\n\t\t"<<"Height of the region from sea level(m)=";
 cin>>h;
 cout<<"\n\t\tThe average Temperature of the region('C)=";
 cin>>T1;
 T1=T1+273;
 //T1 is the average temp of the region for the whole year from previous record
 double P1=Patm1-(h*d*g);
 //P1 is the normal pressure at the region 
 cout<<"\n\t\tPressure observed in the barometer due to weather change is(mm of hg)=";
 double Patm2;
 //Patm2 is the sudden change in pressure at sea level
 cin>>Patm2;
 Patm2=(Patm2*101325)/760;
 double T2;
 //T2 is the Temperature of the region due to pressure change at sea level
 double P2=Patm2-(h*d*g);
 //T2 is the Temperature of the region due to the pressure change at sea level
 T2=T1*(P2/P1);
 if((T2-T1)>0)
 {
     cout<<"\n\t\tThe Temperature('C) of the region has been increased to "<<(T2-273)<<endl;
 }
 else
 {
     cout<<"\n\t\tThe Temperature('C) of the region has been decreased to "<<(T2-273)<<endl;
 }
 cout<<endl<<"\t\tTo continue for further computations press 1"<<endl<<"\t\t";
 cin>>x;
 }while(x==1);
} 
//This is mainly used to convert Temperatures scales to other scales
//Well known scales are Celsius,Fahrenheit,kelvin
//Rankine,Reaumur are industrail notations[not used much]
int MODULE7()                                                                                                              
{                                                                                                                
   temp_conversions();
}
int temp_conversions()

{
    int x;
    do
    { 
        std::system("clear");  
    int n;
    cout<<"\n\t\t\tTEMPERATURE INTERCONVERSIONS";

    cout<<"\n\n\t\t\t1.Interconversion between Celsius and Fahrenheit";
    cout<<"\n\t\t\t2.Interconversion between Celsius and Kelvin";
    cout<<"\n\t\t\t3.Interconversion between fahrenheit and kelvin";
    cout<<"\n\t\t\t4.Design a New Scale for Measuring Temperature";
    cout<<"\n\n\t\tEnter your choice:";
    cin>>n;
    switch(n)
    {
        case 1:
        cel_fahren();
        break;
        case 2:
        cel_kelvin();
        break;
        case 3:
        fahren_kelvin();
        break;
        case 4:
        converting_into_newscale();
        break;
        default:
        cout<<"\n\t\tChoose the right option from 1-4"<<endl;
        break;
    }
    cout<<"\n\t\tTo go to main page press 1:";
    cin>>x;
    }while(x==1);
}
int cel_fahren()
{
    int x;
    do
    {
    float c,f;
    int n;
    cout<<"\n\t\t\t1.Convert from Celsius to Fahrenheit\n\t\t\t2.Convert from Fahrenheit to Celsius";
    cout<<"\n\t\tEnter your choice:";
    cin>>n;
    if(n==1)
    {
        cout<<"\n\n\t\tCELSIUS TO FAHRENHEIT";
        cout<<"\n\tEnter the Temp in Celsius scale=";
        cin>>c;
        f=((9*c)/5)+32;
        cout<<"\n\tTemp in Fahrenheit scale is  "<<f<<endl;
    }
    else
    {
        cout<<"\n\n\t\tFAHRENHEIT TO CELSIUS";
        cout<<"\n\tEnter the Temp in Fahrenheit scale=";
        cin>>f;
        c=(5*(f-32))/9;
        cout<<"\n\tTemp in Celsius scale is  "<<c<<endl;
    }
    cout<<"\n\n\tTo Repeat the computaion press 1:";
    cin>>x;
  }while(x==1);
}
int cel_kelvin()
{
    int x;
    do
    {
    float c,k;
    int n;
    cout<<"\n\t\t\t1.Convert from celsius to kelvin\n\t\t\t2.Convert from Kelvin to Celsius";
    cout<<"\n\t\tEnter your choice:";
    cin>>n;
    if(n==1)
    {
        cout<<"\n\n\t\tCELSIUS TO KELVIN";
        cout<<"\n\tEnter the Temp in Celsius scale=";
        cin>>c;
        k=c+273;
        cout<<"\n\tTemp in Kelvin scale is  "<<k<<endl;
    }
    else
    {
        cout<<"\n\t\t\tKELVIN TO CELSIUS";
        cout<<"\n\tEnter the Temp in Kelvin scale=";
        cin>>k;
        c=k-273;
        cout<<"\n\tTemp in Celsius scale is  "<<c<<endl;
    }
    cout<<"\n\t\tTo Repeat the computation press 1:";
    cin>>x;
   }while(x==1);
}
int fahren_kelvin()
{
    int x;
    do
    {
    float f,k;
    int n;
    cout<<"\n\t\t\t1.Convert from Fahrenheit to Kelvin\n\t\t\t2.Convert from Kelvin to Fahrenheit";
    cout<<"\n\t\tEnter your choice:";
    cin>>n;
    if(n==1)
    {
        cout<<"\n\n\t\tFAHRENHEIT TO KELVIN";
        cout<<"\n\tEnter the Temp in Fahrenheit scale=";
        cin>>f;
        k=((5*(f-32))/9)+273;
        cout<<"\n\tTemp in Kelvin scale is  "<<k<<endl;
    }
    else
    {
        cout<<"\n\n\t\tKELVIN TO FAHRENHEIT";
        cout<<"\n\tEnter the Temp in Kelvin scale=";
        cin>>k;
        f=((9*(k-273))/5)+32;
        cout<<"\n\tTemp in Fahrenheit scale is  "<<f<<endl;
    }
    cout<<"\n\t\tTo Repeat the computation press 1:";
    cin>>x;
   }while(x==1);
}
//This following function is basically used
//to create a new scale for our convenience
//and let's us know waht would be the 
//the value of known scale reading 
//in our new scale
int converting_into_newscale()
{
    int x;
    do
    {   
        cout<<"\n\t\t\t\tINTERCONVERSION FROM KNOWN SCALE TO NEW SCALE";
        int l,u,n,ans,data;
        cout<<"\n\n\t\tEnter the lower and Upper readings of New scale you want to design";
        cout<<"\n\t\tLower reading:";
        cin>>l;
        cout<<"\n\tUpper reading:";
        cin>>u;
        cout<<"\n\n\t1.Convet Celsius to New scale"<<endl;
        cout<<"\n\t2.Convert Fahrenheit to New scale"<<endl;
        cout<<"\n\t3.Convert Kelvin to New scale"<<endl;
        cout<<"\n\tEnter the choice:";
        cin>>n;
        switch(n)
        {
         case 1:
         {
          cout<<"\n\n\t\t1.Convet Celsius to New scale";   
          cout<<"\n\tEnter the Temperature in Celsius scale('C)=";
          cin>>data;
          ans=((data*(u+l))/100)+l;
          cout<<"\n\tThe Temperature in New scale corresponding to Temperature in Celsius is" <<ans<<endl;
          break;
         }
         case 2:
         {
          cout<<"\n\t\t\t2.Convert Fahrenheit to New scale"<<endl;   
          cout<<"\n\tEnter The Temperature in Fahrenheit scale=";
          cin>>data;
          ans=(((data-32)*(u-l))/180)+l;
          cout<<"\n\tThe Temperature in New scale corresponding to Temperature in Fahrenheit is "<<ans<<endl;
          break;
         }
         case 3:
         {
         cout<<"\n\t\t\t3.Convert Kelvin to New scale"<<endl; 
         cout<<"\n\tEnter the Temperature in Kelvin scale(k)=";
         cin>>data;
         ans=(((data-873)*(u-l))/180)+l;
         cout<<"\n\tThe Temperature in New scale corresponding to Temp in Kelvin is "<<ans<<endl;
         break;
         }
         default:
         {
             cout<<"\n\t\tEnter valid number from 1-3"<<endl;
         }
        }
       cout<<"\n\n\t\tIf you want to continue the operation press 1:";
       cin>>x;
    }while(x==1);
}



