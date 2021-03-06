/** ###################################################################
**     THIS BEAN MODULE IS GENERATED BY THE TOOL. DO NOT MODIFY IT.
**     Filename  : Bit1.H
**     Project   : BalanceCar
**     Processor : MCF52255CAF80
**     Beantype  : BitIO
**     Version   : Bean 02.073, Driver 01.01, CPU db: 3.00.000
**     Compiler  : CodeWarrior MCF C Compiler
**     Date/Time : 2015-10-23, 20:03
**     Abstract  :
**         This bean "BitIO" implements an one-bit input/output.
**         It uses one bit/pin of a port.
**         Note: This bean is set to work in Output direction only.
**         Methods of this bean are mostly implemented as a macros
**         (if supported by target language and compiler).
**     Settings  :
**         Used pin                    :
**             ----------------------------------------------------
**                Number (on package)  |    Name
**             ----------------------------------------------------
**                       75            |  FEC_CRS_PTI1
**             ----------------------------------------------------
**
**         Port name                   : PORTTI
**
**         Bit number (in port)        : 1
**         Bit mask of the port        : 0x0002
**
**         Initial direction           : Output (direction cannot be changed)
**         Initial output value        : 1
**         Initial pull option         : up
**
**         Port data register          : PORTTI    [0x00100004]
**         Port control register       : DDRTI     [0x0010001C]
**         Port function register      : PTIPAR    [0x00100064]
**
**         Optimization for            : speed
**     Contents  :
**         GetVal - bool Bit1_GetVal(void);
**         PutVal - void Bit1_PutVal(bool Val);
**         ClrVal - void Bit1_ClrVal(void);
**         SetVal - void Bit1_SetVal(void);
**
**     Copyright : 1997 - 2009 Freescale Semiconductor, Inc. All Rights Reserved.
**     
**     http      : www.freescale.com
**     mail      : support@freescale.com
** ###################################################################*/

#ifndef Bit1_H_
#define Bit1_H_

/* MODULE Bit1. */

  /* Including shared modules, which are used in the whole project */
#include "PE_Types.h"
#include "PE_Error.h"
#include "PE_Const.h"
#include "IO_Map.h"
#include "Cpu.h"


/*
** ===================================================================
**     Method      :  Bit1_GetVal (bean BitIO)
**
**     Description :
**         This method returns an input value.
**           a) direction = Input  : reads the input value from the
**                                   pin and returns it
**           b) direction = Output : returns the last written value
**         Note: This bean is set to work in Output direction only.
**     Parameters  : None
**     Returns     :
**         ---             - Input value. Possible values:
**                           FALSE - logical "0" (Low level)
**                           TRUE - logical "1" (High level)

** ===================================================================
*/
#define Bit1_GetVal() ( \
    (bool)((getReg8(PORTTI) & 0x02))   /* Return port output data */ \
  )

/*
** ===================================================================
**     Method      :  Bit1_PutVal (bean BitIO)
**
**     Description :
**         This method writes the new output value.
**     Parameters  :
**         NAME       - DESCRIPTION
**         Val             - Output value. Possible values:
**                           FALSE - logical "0" (Low level)
**                           TRUE - logical "1" (High level)
**     Returns     : Nothing
** ===================================================================
*/
void Bit1_PutVal(bool Val);

/*
** ===================================================================
**     Method      :  Bit1_ClrVal (bean BitIO)
**
**     Description :
**         This method clears (sets to zero) the output value.
**     Parameters  : None
**     Returns     : Nothing
** ===================================================================
*/
#define Bit1_ClrVal() ( \
    (void)setReg8(CLRTI, 0xFD)         /* CLRTI1=0x00 */ \
  )

/*
** ===================================================================
**     Method      :  Bit1_SetVal (bean BitIO)
**
**     Description :
**         This method sets (sets to one) the output value.
**     Parameters  : None
**     Returns     : Nothing
** ===================================================================
*/
#define Bit1_SetVal() ( \
    (void)setReg8(SETTI, 0x02)         /* SETTI1=0x01 */ \
  )



/* END Bit1. */
#endif /* #ifndef __Bit1_H_ */
/*
** ###################################################################
**
**     This file was created by Processor Expert 1.05 [04.27]
**     for the Freescale MCF series of microcontrollers.
**
** ###################################################################
*/
