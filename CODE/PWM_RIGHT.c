/** ###################################################################
**     THIS BEAN MODULE IS GENERATED BY THE TOOL. DO NOT MODIFY IT.
**     Filename  : PWM_RIGHT.C
**     Project   : BalanceCar
**     Processor : MCF52255CAF80
**     Beantype  : PWM
**     Version   : Bean 02.210, Driver 01.04, CPU db: 3.00.000
**     Compiler  : CodeWarrior MCF C Compiler
**     Date/Time : 2014-12-6, 21:42
**     Abstract  :
**         This bean implements a pulse-width modulation generator
**         that generates signal with variable duty and fixed cycle. 
**     Settings  :
**         Used output pin             : 
**             ----------------------------------------------------
**                Number (on package)  |    Name
**             ----------------------------------------------------
**                       5             |  DTIN0_DTOUT0_PWM0_PTC0
**             ----------------------------------------------------
**
**         Timer name                  : PWMCNT0 [8-bit]
**         Counter                     : PWMCNT0   [0x001B000C]
**         Mode register               : PWMCTL    [0x001B0005]
**         Run register                : PWME      [0x001B0000]
**         Prescaler                   : PWMPRCLK  [0x001B0003]
**         Compare 1 register          : PWMPER0   [0x001B0014]
**         Compare 2 register          : PWMDTY0   [0x001B001C]
**         Flip-flop 1 register        : PWMPOL    [0x001B0001]
**
**         User handling procedure     : not specified
**
**         Output pin
**
**         Port name                   : PORTTC
**         Bit number (in port)        : 0
**         Bit mask of the port        : 0x0001
**         Port data register          : PORTTC    [0x0010000F]
**         Port control register       : DDRTC     [0x00100027]
**
**         Runtime setting period      : none
**         Runtime setting ratio       : calculated
**         Initialization:
**              Aligned                : Left
**              Output level           : low
**              Timer                  : Enabled
**              Event                  : Enabled
**         High speed mode
**             Prescaler               : divide-by-1
**             Clock                   : 2500000 Hz
**           Initial value of            period        pulse width (ratio 50%)
**             Xtal ticks              : 8000          4000
**             microseconds            : 100           50
**             seconds (real)          : 0.0001        0.00005
**
**     Contents  :
**         SetRatio16 - byte PWM_RIGHT_SetRatio16(word Ratio);
**         SetDutyUS  - byte PWM_RIGHT_SetDutyUS(word Time);
**         SetDutyMS  - byte PWM_RIGHT_SetDutyMS(word Time);
**
**     Copyright : 1997 - 2009 Freescale Semiconductor, Inc. All Rights Reserved.
**     
**     http      : www.freescale.com
**     mail      : support@freescale.com
** ###################################################################*/

/* MODULE PWM_RIGHT. */

#include "PWM_RIGHT.h"

static word RatioStore;                /* Ratio of L-level to H-level */
/*
** ===================================================================
**     Method      :  SetRatio (bean PWM)
**
**     Description :
**         The method reconfigures the compare and modulo registers of 
**         the peripheral(s) when the speed mode changes. The method is 
**         called automatically as a part of the bean 
**         SetHigh/SetLow/SetSlow methods.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/
static void SetRatio(void)
{
  /* Calculate new value according to the given ratio */
  setReg8(PWMDTY0,(byte)(((((dword)PWMPER0 + 1) * RatioStore) + (PWMPER0 & 0x8000)) >> 0x10));
  /* PWMCNT0: COUNT=0 */
  setReg8(PWMCNT0, 0x00);              /* Reset counter */ 
}

/*
** ===================================================================
**     Method      :  PWM_RIGHT_SetRatio16 (bean PWM)
**
**     Description :
**         This method sets a new duty-cycle ratio.
**     Parameters  :
**         NAME       - DESCRIPTION
**         Ratio      - Ratio is expressed as an 16-bit unsigned integer
**                      number. 0 - 0xFFFF value is proportional
**                      to ratio 0 - 100%
**         Note: Calculated duty depends on the timer possibilities
**               and on the selected period.
**     Returns     :
**         ---        - Error code, possible codes:
**                           ERR_OK - OK
**                           ERR_SPEED - This device does not work in
**                           the active speed mode
** ===================================================================
*/
byte PWM_RIGHT_SetRatio16(word Ratio)
{
  RatioStore = Ratio;                  /* Store new value of the ratio */
  SetRatio();                          /* Calculate and set up new appropriate values of the duty and period registers */
  return ERR_OK;                       /* OK */
}

/*
** ===================================================================
**     Method      :  PWM_RIGHT_SetDutyUS (bean PWM)
**
**     Description :
**         This method sets the new duty value of the output signal. The
**         duty is expressed in microseconds as a 16-bit unsigned integer
**         number.
**     Parameters  :
**         NAME       - DESCRIPTION
**         Time       - Duty to set [in microseconds]
**                      (0 to 100 us in high speed mode)
**     Returns     :
**         ---        - Error code, possible codes:
**                           ERR_OK - OK
**                           ERR_SPEED - This device does not work in
**                           the active speed mode
**                           ERR_MATH - Overflow during evaluation
**                           ERR_RANGE - Parameter out of range
** ===================================================================
*/
byte PWM_RIGHT_SetDutyUS(word Time)
{
  dlong rtmul;                         /* Result of two 32-bit numbers multiplication */
  
  if (Time > 0x64) {                   /* Is the given value out of range? */
    return ERR_RANGE;                  /* If yes then error */
  }
  rtmul = (dlong)Time * 0x028F5C29UL;  /* Multiply given value and High speed CPU mode coefficient */
  if (PE_Timer_LngHi2(rtmul,&RatioStore)) { /* Is the result greater or equal than 65536 ? */
    RatioStore = 0xFFFF;               /* If yes then use maximal possible value */
  }
  SetRatio();                          /* Calculate and set up new appropriate values of the duty and period registers */
  return ERR_OK;                       /* OK */
}

/*
** ===================================================================
**     Method      :  PWM_RIGHT_SetDutyMS (bean PWM)
**
**     Description :
**         This method sets the new duty value of the output signal. The
**         duty is expressed in milliseconds as a 16-bit unsigned integer
**         number.
**     Parameters  :
**         NAME       - DESCRIPTION
**         Time       - Duty to set [in milliseconds]
**         Note: The period is too short. The method will return
**               just the error code in high speed mode.
**     Returns     :
**         ---        - Error code, possible codes:
**                           ERR_OK - OK
**                           ERR_SPEED - This device does not work in
**                           the active speed mode
**                           ERR_MATH - Overflow during evaluation
**                           ERR_RANGE - Parameter out of range
** ===================================================================
*/
#pragma warn_unusedarg off
byte PWM_RIGHT_SetDutyMS(word Time)
{
  /* Period is too little. Method 'SetDutyMS' will return only error code in High speed CPU mode. */
  return ERR_MATH;                     /* Calculation error */
}
#pragma warn_unusedarg reset

/*
** ===================================================================
**     Method      :  PWM_RIGHT_Init (bean PWM)
**
**     Description :
**         Initializes the associated peripheral(s) and the beans 
**         internal variables. The method is called automatically as a 
**         part of the application initialization code.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/
void PWM_RIGHT_Init(void)
{
  /* PWMDTY0: DUTY=0x7D */
  setReg16(PWMDTY0, 0x7DU);            /* Store initial value to the duty-compare register */ 
  /* PWMPER0: PERIOD=0xFA */
  setReg16(PWMPER0, 0xFAU);            /* and to the period register */ 
  /* PWMCNT0: COUNT=0 */
  setReg16(PWMCNT0, 0x00U);            /* Reset counter */ 
  RatioStore = 0x7F7E;
}

/* END PWM_RIGHT. */

/*
** ###################################################################
**
**     This file was created by Processor Expert 1.05 [04.27]
**     for the Freescale MCF series of microcontrollers.
**
** ###################################################################
*/