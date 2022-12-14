//+------------------------------------------------------------------+
//|                                       GeoMotionDemonstration.mq4 |
//|                                                      WJBOSS INC. |
//|                                                             None |
//+------------------------------------------------------------------+
#property strict
#property indicator_minimum -20
#property indicator_maximum  20

#property indicator_separate_window

#property indicator_buffers 1

#property indicator_color1 Red
#property indicator_width1  1


input int InpN_Number_Of_Computations =1000;
double ExtSignalBuffer[];
bool ExtParameters =false;

input ENUM_TIMEFRAMES TimeFrame = PERIOD_M1;

input int NumberOfStpes = 120;

input int CandlesIntoTheFuture= 60; 


//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
  {
  IndicatorDigits(Digits+1);
  
  
  SetIndexStyle(0,DRAW_LINE);
  
  SetIndexBuffer(0,ExtSignalBuffer);
  
  IndicatorShortName("GEO:("+IntegerToString(InpN_Number_Of_Computations)+")");
  
  SetIndexLabel(0,"GEO");
  
  if(InpN_Number_Of_Computations < 100){
  
  Print("too few computations");
  ExtParameters = false;
  return(INIT_FAILED);
  
  }
  else{
  
  ExtParameters = true;
  return(INIT_SUCCEEDED);
  
  
  }
   
//---
   return(INIT_SUCCEEDED);
  }
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &time[],
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
  {
//---
   int i, limit;
   
   limit = rates_total-prev_calculated;
   
   if( prev_calculated>0){
   
   limit++;}
   
   
   
   int Size = ArraySize(ExtSignalBuffer);
   
   ArrayInitialize(ExtSignalBuffer,0);
   
   for(i = 0 ;i<61; i++){
   
   
   
      double mu = Point();
      
      int n = NumberOfStpes;
      double T = CandlesIntoTheFuture*double(TimeFrame)/525600.0;
      int T_Hat = CandlesIntoTheFuture*double(TimeFrame);
      
      
      int M = InpN_Number_Of_Computations;
      
      double S0 = close[T_Hat-1 +i];
      
      double Sigma = 0;
      
      
      for( int j = T_Hat-1 ; j>=1; j--){
      
      
      Sigma += pow(  MathLog(close[i+j-1]/close[i+j]),2);
      
      
      
      }
      
      Sigma = 100*MathSqrt((1/T_Hat)*Sigma);
       
      double dt = T/double(n);
      
      double Geo = GeometericBrownianMotion(S0,TimeFrame,mu,Sigma,NumberOfStpes,dt,InpN_Number_Of_Computations);
      
      ExtSignalBuffer[i] = (Geo) - close[i]; 
   
   
   
   
   }
//--- return value of prev_calculated for next call
   return(rates_total);
  }
//+------------------------------------------------------------------+


double rand_gen(){




   return ( (double)(rand())+1.0)/((double)(32767.0)+1.0);

}


double normalRandom(){


   double v1 = rand_gen();
   double v2 = rand_gen();
   
   return  cos(2*3.141592653589793*v2)*sqrt(-2.0*log(v1));


}


double NormalDistSample(double Mu, double sigma){



   double x = normalRandom()*sigma+Mu;
   return x;
   



}


double GeometericBrownianMotion(double S0,ENUM_TIMEFRAMES TF,double mu, double sigma,int steps ,double dt,int n_computations){



double ST_0 = 0;

double ST_mu = 0;


   for( int j = 0 ; j < n_computations; j++){


      ST_0= S0;
      
      for( int i =0; i < steps; i++){
         
         ST_0 = ST_0*MathExp(  (mu - .5*pow(sigma,2))*dt+sigma*(NormalDistSample(0,MathSqrt(dt)))); 
      
      }


      ST_mu += ST_0;

   }

   ST_mu = ST_mu/double(n_computations);
   return ST_mu;

}