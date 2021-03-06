/** ###################################################################
**     THIS BEAN MODULE IS GENERATED BY THE TOOL. DO NOT MODIFY IT.
**     Filename  : PE_Const.H
**     Project   : BalanceCar
**     Processor : MCF52255CAF80
**     Beantype  : PE_Const
**     Version   : Driver 01.00
**     Compiler  : CodeWarrior MCF C Compiler
**     Date/Time : 2013/9/28, 21:08
**     Abstract  :
**         This bean "PE_Const" contains internal definitions
**         of the constants.
**     Settings  :
**     Contents  :
**         No public methods
**
**     Copyright : 1997 - 2009 Freescale Semiconductor, Inc. All Rights Reserved.
**     
**     http      : www.freescale.com
**     mail      : support@freescale.com
** ###################################################################*/

#ifndef __PE_Const_H
#define __PE_Const_H

/* Constants for detecting running mode */
#define HIGH_SPEED        0x00         /* High speed */
#define LOW_SPEED         0x01         /* Low speed */
#define SLOW_SPEED        0x02         /* Slow speed */

/* Reset cause constants */
#define RSTSRC_LOL        0x01         /* Loss-of-lock reset */
#define RSTSRC_LOC        0x02         /* Loss-of-clock reset */
#define RSTSRC_EXT        0x04         /* External reset */
#define RSTSRC_POR        0x08         /* Power-on reset */
#define RSTSRC_SOFT       0x20         /* Software reset */
#define RSTSRC_LVD        0x40         /* Low voltage detect reset */

#endif /* _PE_Const_H */
/*
** ###################################################################
**
**     This file was created by Processor Expert 1.05 [04.27]
**     for the Freescale MCF series of microcontrollers.
**
** ###################################################################
*/
