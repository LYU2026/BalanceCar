/** ###################################################################
**     THIS BEAN MODULE IS GENERATED BY THE TOOL. DO NOT MODIFY IT.
**     Filename  : Cpu.C
**     Project   : BalanceCar
**     Processor : MCF52255CAF80
**     Beantype  : MCF52255_100_LQFP
**     Version   : Bean 01.000, Driver 01.02, CPU db: 3.00.000
**     Datasheet : MCF52259RM Rev.0 Draft B, MCF52259DS Rev.0 Draft C, Kirin3 System on a Chip GuideV0.2
**     Compiler  : CodeWarrior MCF C Compiler
**     Date/Time : 2014-11-23, 15:52
**     Abstract  :
**
**     Settings  :
**
**     Contents  :
**         EnableInt   - void Cpu_EnableInt(void);
**         DisableInt  - void Cpu_DisableInt(void);
**         SetIntLevel - void Cpu_SetIntLevel(byte Level);
**         Delay100US  - void Cpu_Delay100US(word us100);
**
**     Copyright : 1997 - 2009 Freescale Semiconductor, Inc. All Rights Reserved.
**     
**     http      : www.freescale.com
**     mail      : support@freescale.com
** ###################################################################*/

/* MODULE Cpu. */
#include "Test.h"
#include "PE_Types.h"
#include "PE_Error.h"
#include "PE_Const.h"
#include "IO_Map.h"
#include "Events.h"
#include "Cpu.h"


/* Global variables */
volatile word SR_reg;                  /* Current value of the SR register */
volatile byte SR_lock = 0;             /* Lock */
/*
** ===================================================================
**     Method      :  Cpu_Interrupt (bean MCF52255_100_LQFP)
**
**     Description :
**         The method services unhandled interrupt vectors.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/
__declspec(interrupt) void Cpu_Interrupt(void)
{
  asm(HALT);
}

/*
** ===================================================================
**     Method      :  Cpu_INT_SCM_SWTIInterrupt (bean MCF52255_100_LQFP)
**
**     Description :
**         This ISR services the 'Core watchdog' interrupt.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/
__declspec(interrupt) void Cpu_INT_SCM_SWTIInterrupt(void)
{
  Cpu_OnCoreWatchdogINT();
}

/*
** ===================================================================
**     Method      :  Cpu_DisableInt (bean MCF52255_100_LQFP)
**
**     Description :
**         Disables all maskable interrupts. This method sets the
**         interrupt level mask group in the SR register by value = 7.
**     Parameters  : None
**     Returns     : Nothing
** ===================================================================
*/
/*
void Cpu_DisableInt(void)

**      This method is implemented as macro in the header module. **
*/

/*
** ===================================================================
**     Method      :  Cpu_EnableInt (bean MCF52255_100_LQFP)
**
**     Description :
**         Enables all maskable interrupts. This method sets the
**         interrupt level mask group in the SR register by value = 0.
**     Parameters  : None
**     Returns     : Nothing
** ===================================================================
*/
/*
void Cpu_EnableInt(void)

**      This method is implemented as macro in the header module. **
*/

/*
** ===================================================================
**     Method      :  Cpu_SetIntLevel (bean MCF52255_100_LQFP)
**
**     Description :
**         Sets the interrupt level mask in the SR register. Interrupt
**         requests are inhibited for all priority levels less than or
**         equal to current level, except edge-sensitive level 7
**         request, which cannot be masked.
**     Parameters  :
**         NAME            - DESCRIPTION
**         Level           - New interrupt priority level value.
**     Returns     : Nothing
** ===================================================================
*/
/*
void Cpu_SetIntLevel(void)

**      This method is implemented as macro in the header module. **
*/

/*
** ===================================================================
**     Method      :  Cpu_Delay100US (bean MCF52255_100_LQFP)
**
**     Description :
**         This method realizes software delay. The length of delay is
**         at least 100 microsecond multiply input parameter [us100].
**         As the delay implementation is not based on real clock, the
**         delay time may be increased by interrupt service routines
**         processed during the delay. The method is independent on
**         selected speed mode.
**     Parameters  :
**         NAME            - DESCRIPTION
**         us100           - Number of 100 us delay repetitions.
**     Returns     : Nothing
** ===================================================================
*/
__declspec(register_abi) void Cpu_Delay100US(word us100:__D0)
{
  /* irremovable one time overhead (ignored): 11 cycles */
  /* move: 1 cycle overhead (load parameter into register) */
  /* jsr:  3 cycles overhead (jump to subroutine) */
  /* andi: 1 cycle overhead (clear upper word of d0) */
  /* tpf: 1 cycle overhead (alignment) */
  /* rts:  5 cycles overhead (return from subroutine) */

  /* aproximate irremovable overhead for each 100us cycle (counted) : 3 cycles */
  /* subq.l 1 cycles overhead  */
  /* bne.b  2 cycles overhead  */

  /* Disable MISRA rule 55 checking - Non-case label used */
  /*lint -esym( 961, 55)   */
#pragma unused(us100)
  asm {
    naked
    andi.l #0xFFFF,d0                  /* parameter is word - clear the rest of d0 register */
    tpf                                /* alignment */
loop:
    /* 100 us delay block begin */
    /*
     * Delay
     *   - requested                  : 100 us @ 80MHz,
     *   - possible                   : 8000 c, 100000 ns
     *   - without removable overhead : 7997 c, 99962.5 ns
     */
    move.l #0x0A69,d1                  /* (1 c: 12.5 ns) number of iterations */
label0:
    subq.l #1,d1                       /* (1 c: 12.5 ns) decrement d1 */
    bne.b label0                       /* (2 c: 25 ns) repeat 2665x */
    tpf                                /* (1 c: 12.5 ns) wait for 1 c */
    /* 100 us delay block end */
    subq.l #1,d0                       /* parameter is passed via d0 register */
    bne.b loop                         /* next loop */
    rts                                /* return from subroutine */
  }
  /* Restore MISRA rule 55 checking - Non-case label used */
  /*lint +esym( 961, 55)   */
}

/*
** ===================================================================
**     Method      :  __initialize_hardware (bean MCF52255_100_LQFP)
**
**     Description :
**         Initializes the whole system like timing, external bus, etc.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/

/*** !!! Here you can place your own code using property "User data declarations" on the build options tab. !!! ***/

void __initialize_hardware(void)
{
  extern uint32 __VECTOR_RAM[];


  /*** !!! Here you can place your own code before PE initialization using property "User code before PE initialization" on the build options tab. !!! ***/


  /* Initialize IPSBAR and validate it */
  asm {
    move.l  #0x40000001,d0
    move.l  d0,0x40000000
  }
  /* Initialize FLASHBAR */
  asm {
    move.l  #0x00000021,d0
    movec   d0,FLASHBAR
  }
  /* Initialize RAMBAR */
  asm {
    move.l  #0x20000221,d0
    movec   d0,RAMBAR
  }
  /* Initialize Vector Base Register (VBR) */
  asm {
    lea   __VECTOR_RAM,A0
    movec A0,VBR
    nop
  }

  /*** ### MCF52255CAF80 "Cpu" init code ... ***/
  /*** PE initialization code after reset ***/
  /* System clock initialization */
  /* CCHR: ??=0,??=0,??=0,??=0,??=0,CCHR=0 */
  setReg8(CCHR, 0x00);                 /* Set the predivider */ 
  /* SYNCR: LOLRE=0,MFD=3,LOCRE=0,RFD=0,LOCEN=0,DISCLK=0,FWKUP=0,??=0,??=0,CLKSRC=1,PLLMODE=1,PLLEN=1 */
  setReg16(SYNCR, 0x3007U);            /* Set the SYNCR register */ 
  while (!(SYNSR & SYNSR_LOCK_BITMASK)){} /* Wait until the PLL is locked. */
  /* LPDR: ??=0,??=0,??=0,??=0,LPD=0 */
  setReg8(LPDR, 0x00);                 /* Set the low power divider */ 
  /* RTCCR: ??=0,EXTALEN=0,??=0,OSCEN=0,KHZEN=0,REFS=0,LPEN=0,RTCSEL=1 */
  setReg8(RTCCR, 0x01);                /* Set the RTC oscillator and RTC clock select */ 
  /* BWCR: ??=0,??=0,??=0,??=0,??=0,??=0,BWDSTOP=0,BWDSEL=1 */
  setReg32(BWCR, 0x01UL);              /* Set the BWCR register */ 
  /* Common initialization of the CPU registers */
  /*** End of PE initialization code after reset ***/

  /*** !!! Here you can place your own code after PE initialization using property "User code after PE initialization" on the build options tab. !!! ***/

}

/*
** ===================================================================
**     Method      :  PE_low_level_init (bean MCF52255_100_LQFP)
**
**     Description :
**         Initializes beans and provides common register initialization. 
**         The method is called automatically as a part of the 
**         application initialization code.
**         This method is internal. It is used by Processor Expert only.
** ===================================================================
*/
void PE_low_level_init(void)
{
  /* Initialization of the SCM module */
  /* CWCR: CWE=0,CWRI=0,CWT=0,CWTA=0,CWTAVAL=1,CWTIF=1 */
  setReg8(CWCR, 0x03);                  
  /* CWSR: CWSR=0x55 */
  setReg8(CWSR, 0x55);                  
  /* CWSR: CWSR=0xAA */
  setReg8(CWSR, 0xAA);                  
  /* CWCR: CWE=0,CWRI=0,CWT=0,CWTA=0,CWTAVAL=0,CWTIF=0 */
  setReg8(CWCR, 0x00);                  
  /* MPARK: BCR24BIT=0 */
  clrReg32Bits(MPARK, 0x01000000UL);    
  /* DMAREQC: DMAC3=0,DMAC2=0,DMAC1=0,DMAC0=0 */
  clrReg32Bits(DMAREQC, 0xFFFFUL);      
  /* MPARK: M2_P_EN=0,M3_PRTY=3,M2_PRTY=2,M0_PRTY=0,M1_PRTY=1,FIXED=0,TIMEOUT=0,PRKLAST=0,LCKOUT_TIME=0 */
  clrSetReg32Bits(MPARK, 0x021E7F00UL, 0x00E10000UL); 
  /* MPR: MPR=3 */
  clrSetReg8Bits(MPR, 0x0C, 0x03);      
  /* PACR0: LOCK1=0,ACCESS_CTRL1=0 */
  clrReg8Bits(PACR0, 0xF0);             
  /* PACR1: LOCK0=0,ACCESS_CTRL0=0 */
  clrReg8Bits(PACR1, 0x0F);             
  /* PACR2: LOCK1=0,ACCESS_CTRL1=0,LOCK0=0,ACCESS_CTRL0=0 */
  setReg8(PACR2, 0x00);                 
  /* PACR3: LOCK1=0,ACCESS_CTRL1=0 */
  clrReg8Bits(PACR3, 0xF0);             
  /* PACR4: LOCK1=0,ACCESS_CTRL1=0,LOCK0=0,ACCESS_CTRL0=0 */
  setReg8(PACR4, 0x00);                 
  /* PACR5: LOCK1=0,ACCESS_CTRL1=0 */
  clrReg8Bits(PACR5, 0xF0);             
  /* PACR6: LOCK1=0,ACCESS_CTRL1=0,LOCK0=0,ACCESS_CTRL0=0 */
  setReg8(PACR6, 0x00);                 
  /* PACR7: LOCK1=0,ACCESS_CTRL1=0,LOCK0=0,ACCESS_CTRL0=0 */
  setReg8(PACR7, 0x00);                 
  /* PACR8: LOCK1=0,ACCESS_CTRL1=0,LOCK0=0,ACCESS_CTRL0=0 */
  setReg8(PACR8, 0x00);                 
  /* PACR10: LOCK1=0,ACCESS_CTRL1=0 */
  clrReg8Bits(PACR10, 0xF0);            
  /* GPACR0: LOCK=0,ACCESS_CTRL=0 */
  clrReg8Bits(GPACR0, 0x8F);            
  /* GPACR1: LOCK=0,ACCESS_CTRL=0 */
  clrReg8Bits(GPACR1, 0x8F);            
  /* Initialization of the PowerManagement module */
  /* LPCR: LPMD=0 */
  clrReg8Bits(LPCR, 0xC0);              
  /* LPICR: ENBSTOP=0,XLPM_IPL=0 */
  clrReg8Bits(LPICR, 0xF0);             
  /* PPMRH: CDUSB=0,CDCFM=0,CDPWM=0,CDGPT=0,CDADC=0,CDPIT1=0,CDPIT0=0,CDEPORT=0,CDGPIO=0 */
  clrReg32Bits(PPMRH, 0x1B9BUL);        
  /* PPMRL: CDFEC=0,CDINTC1=0,CDINTC0=0,CDDTIM3=0,CDDTIM2=0,CDDTIM1=0,CDDTIM0=0,CDRTC=0,CDI2C1=0,CDQSPI=0,CDI2C0=0,CDUART2=0,CDUART1=0,CDUART0=0,CDDMA=0,CDG=0 */
  clrReg32Bits(PPMRL, 0x0027FEF2UL);    
  /* IPSBMT: BME=1,BMT=0 */
  clrSetReg8Bits(IPSBMT, 0x07, 0x08);   
  /* Initialization of the ResetController module */
  /* ResetController_RCR: SOFTRST=0,FRCRSTOUT=0,??=0,LVDF=1,LVDIE=0,LVDRE=1,??=0,LVDE=1 */
  setReg8(ResetController_RCR, 0x15);   
  /* Initialization of the CCM module */
  /* CCR: LOAD=0 */
  clrReg16Bits(CCR, 0x8000U);           
  /* IMRL0: MASKALL=0 */
  clrReg32Bits(IMRL0, 0x01UL);         /* Disable masking of all interrupts of the INTC0 device */ 
  /* IMRL1: MASKALL=0 */
  clrReg32Bits(IMRL1, 0x01UL);         /* Disable masking of all interrupts of the INTC1 device */ 
  /* SCM_RAMBAR: BA=0x2000,??=0,??=0,??=0,??=0,??=0,??=0,BDE=1,??=0,??=0,??=0,??=0,??=0,??=0,??=0,??=0,??=0 */
  setReg32(SCM_RAMBAR, 0x20000200UL);   
  /* Common initialization of the CPU registers */
  /* PORTTI: PORTTI0=1 */
  setReg8Bits(PORTTI, 0x01);            
  /* DDRTI: DDRTI0=1 */
  setReg8Bits(DDRTI, 0x01);             
  /* PTIPAR: PTIPAR0=0 */
  clrReg8Bits(PTIPAR, 0x01);            
  /* ICR008: IL=0,IP=0 */
  clrReg8Bits(ICR008, 0x3F);            
  /* ICR048: IL=0,IP=0 */
  clrReg8Bits(ICR048, 0x3F);            
  /* IMRL0: INT_MASK8=1 */
  setReg32Bits(IMRL0, 0x0100UL);        
  /* IMRH0: INT_MASK48=1 */
  setReg32Bits(IMRH0, 0x00010000UL);    
  /* ### BitIO "Test" init code ... */
  __EI(0);                             /* Enable interrupts of the given priority level */
}

/* END Cpu. */

/*
** ###################################################################
**
**     This file was created by Processor Expert 1.05 [04.27]
**     for the Freescale MCF series of microcontrollers.
**
** ###################################################################
*/
