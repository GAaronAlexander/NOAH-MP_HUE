  LY  ^   k820309    ?          19.1        Ïe-e                                                                                                          
       module_sf_noahmp_groundwater.f90 MODULE_SF_NOAHMP_GROUNDWATER #         @                                                    0   #NSOIL    #XLAND    #XICE    #XICE_THRESHOLD 	   #ISICE 
   #ISLTYP    #SMOISEQ    #DZS    #WTDDT    #FDEPTH    #AREA    #TOPO    #ISURBAN    #IVGTYP    #RIVERCOND    #RIVERBED    #EQWTD    #PEXP    #SMOIS    #SH2OXY    #SMCWTD    #WTD    #QLAT    #QRF    #DEEPRECH    #QSPRING    #QSLAT     #QRFS !   #QSPRINGS "   #RECH #   #IDS $   #IDE %   #JDS &   #JDE '   #KDS (   #KDE )   #IMS    #IME    #JMS    #JME    #KMS *   #KME +   #ITS ,   #ITE -   #JTS .   #JTE /   #KTS 0   #KTE 1             
  @                                                   
                                                      	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                    
                                  	     	                
                                  
                    
  @                                                         5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	        5  p '       r    5  p        r    p          5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & p        5  p        r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p          5  p        r        5  p (       r    5  p '       r    p                                   
  @                                                   	    p          & p        5  p        r        5  p        r                                
                                       	               
  @                                                   	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
  @                                                   	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
  @                                                   	 	     5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                    
                                                      
                                                            5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	 
     5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
                                                      	      5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                                     	         5  p '       r    5  p        r    p          5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & p        5  p        r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p          5  p        r        5  p (       r    5  p '       r    p                                   
D                                                     	         5  p '       r    5  p        r    p          5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & p        5  p        r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p          5  p        r        5  p (       r    5  p '       r    p                                   
D @                                                   	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D @                                                   	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   D @                                                   	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   D                                                     	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                                     	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   D @                                                   	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                                      	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                 !                    	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                 "                    	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                   
D                                 #                    	       5  p '       r      5  p &       r    5  p %       r    p        5  p %       r      & 5  p %       r    5  p &       r      & 5  p '       r    5  p (       r          5  p &       r    5  p %       r    p            5  p (       r    5  p '       r    p                                    
  @                               $                     
  @                               %                     
  @                               &                     
  @                               '                     
  @                               (                     
  @                               )                     
  @                                                    
  @                                                    
  @                                                    
  @                                                    
  @                               *                     
  @                               +                     
  @                               ,                     
  @                               -                     
  @                               .                     
  @                               /                     
  @                               0                     
  @                               1           #         @                                  2                    #ISLTYP 3   #WTD 8   #QLAT 9   #FDEPTH :   #TOPO ;   #LANDMASK <   #DELTAT =   #AREA >   #IDS ?   #IDE @   #JDS A   #JDE B   #KDS C   #KDE D   #IMS 6   #IME 5   #JMS 4   #JME 7   #KMS E   #KME F   #ITS G   #ITE H   #JTS I   #JTE J   #KTS K   #KTE L            
                                  3                          5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                   
                                  8                    	 "     5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                   D                                 9                    	 %      5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                   
                                  :                    	 !     5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                   
                                  ;                    	 #     5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                   
                                  <                           5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                    
                                  =     	               
                                  >                    	 $     5  p        r 4     5  p        r 5   5  p        r 6   p        5  p        r 6     & 5  p        r 6   5  p        r 5     & 5  p        r 4   5  p        r 7         5  p        r 5   5  p        r 6   p            5  p        r 7   5  p        r 4   p                                    
  @                               ?                     
                                  @                     
  @                               A                     
                                  B                     
                                  C                     
                                  D                     
                                  6                     
                                  5                     
                                  4                     
                                  7                     
                                  E                     
                                  F                     
  @                               G                     
  @                               H                     
  @                               I                     
  @                               J                     
                                  K                     
                                  L           #         @                                  M                    #NSOIL N   #DZS O   #ZSOIL P   #SMCEQ Q   #SMCMAX R   #SMCWLT S   #PSISAT T   #BEXP U   #ILOC V   #JLOC W   #TOTWATER X   #WTD Y   #SMC Z   #SH2O [   #SMCWTD \   #QSPRING ]             
                                  N                    
                                  O                    	 +   p          & p        5  p        r N       5  p        r N                              
                                  P                    	 )   p           & p         5  p        r N         5  p        r N   p         p                                   
                                  Q                    	 *   p          & p        5  p        r N       5  p        r N                               
  @                               R     	                
                                  S     	                
                                  T     	                
                                  U     	                
                                  V                     
                                  W                     
D                                 X     	                 
D                                 Y     	                
D                                 Z                    	 ,    p          & p        5  p        r N       5  p        r N                              
D                                 [                    	 -    p          & p        5  p        r N       5  p        r N                               
D @                               \     	                 D                                 ]     	              F      fn#fn "   æ   H      WTABLE_MMF_NOAHMP (   .  @   a   WTABLE_MMF_NOAHMP%NSOIL (   n    a   WTABLE_MMF_NOAHMP%XLAND '       a   WTABLE_MMF_NOAHMP%XICE 1     @   a   WTABLE_MMF_NOAHMP%XICE_THRESHOLD (   Ö  @   a   WTABLE_MMF_NOAHMP%ISICE )       a   WTABLE_MMF_NOAHMP%ISLTYP *   *
    a   WTABLE_MMF_NOAHMP%SMOISEQ &   ¾  Ä   a   WTABLE_MMF_NOAHMP%DZS (     @   a   WTABLE_MMF_NOAHMP%WTDDT )   Â    a   WTABLE_MMF_NOAHMP%FDEPTH '   Ö    a   WTABLE_MMF_NOAHMP%AREA '   ê    a   WTABLE_MMF_NOAHMP%TOPO *   þ  @   a   WTABLE_MMF_NOAHMP%ISURBAN )   >    a   WTABLE_MMF_NOAHMP%IVGTYP ,   R    a   WTABLE_MMF_NOAHMP%RIVERCOND +   f    a   WTABLE_MMF_NOAHMP%RIVERBED (   z    a   WTABLE_MMF_NOAHMP%EQWTD '       a   WTABLE_MMF_NOAHMP%PEXP (   ¢    a   WTABLE_MMF_NOAHMP%SMOIS )   6!    a   WTABLE_MMF_NOAHMP%SH2OXY )   Ê#    a   WTABLE_MMF_NOAHMP%SMCWTD &   Þ%    a   WTABLE_MMF_NOAHMP%WTD '   ò'    a   WTABLE_MMF_NOAHMP%QLAT &   *    a   WTABLE_MMF_NOAHMP%QRF +   ,    a   WTABLE_MMF_NOAHMP%DEEPRECH *   ..    a   WTABLE_MMF_NOAHMP%QSPRING (   B0    a   WTABLE_MMF_NOAHMP%QSLAT '   V2    a   WTABLE_MMF_NOAHMP%QRFS +   j4    a   WTABLE_MMF_NOAHMP%QSPRINGS '   ~6    a   WTABLE_MMF_NOAHMP%RECH &   8  @   a   WTABLE_MMF_NOAHMP%IDS &   Ò8  @   a   WTABLE_MMF_NOAHMP%IDE &   9  @   a   WTABLE_MMF_NOAHMP%JDS &   R9  @   a   WTABLE_MMF_NOAHMP%JDE &   9  @   a   WTABLE_MMF_NOAHMP%KDS &   Ò9  @   a   WTABLE_MMF_NOAHMP%KDE &   :  @   a   WTABLE_MMF_NOAHMP%IMS &   R:  @   a   WTABLE_MMF_NOAHMP%IME &   :  @   a   WTABLE_MMF_NOAHMP%JMS &   Ò:  @   a   WTABLE_MMF_NOAHMP%JME &   ;  @   a   WTABLE_MMF_NOAHMP%KMS &   R;  @   a   WTABLE_MMF_NOAHMP%KME &   ;  @   a   WTABLE_MMF_NOAHMP%ITS &   Ò;  @   a   WTABLE_MMF_NOAHMP%ITE &   <  @   a   WTABLE_MMF_NOAHMP%JTS &   R<  @   a   WTABLE_MMF_NOAHMP%JTE &   <  @   a   WTABLE_MMF_NOAHMP%KTS &   Ò<  @   a   WTABLE_MMF_NOAHMP%KTE    =  C      LATERALFLOW #   U>    a   LATERALFLOW%ISLTYP     i@    a   LATERALFLOW%WTD !   }B    a   LATERALFLOW%QLAT #   D    a   LATERALFLOW%FDEPTH !   ¥F    a   LATERALFLOW%TOPO %   ¹H    a   LATERALFLOW%LANDMASK #   ÍJ  @   a   LATERALFLOW%DELTAT !   K    a   LATERALFLOW%AREA     !M  @   a   LATERALFLOW%IDS     aM  @   a   LATERALFLOW%IDE     ¡M  @   a   LATERALFLOW%JDS     áM  @   a   LATERALFLOW%JDE     !N  @   a   LATERALFLOW%KDS     aN  @   a   LATERALFLOW%KDE     ¡N  @   a   LATERALFLOW%IMS     áN  @   a   LATERALFLOW%IME     !O  @   a   LATERALFLOW%JMS     aO  @   a   LATERALFLOW%JME     ¡O  @   a   LATERALFLOW%KMS     áO  @   a   LATERALFLOW%KME     !P  @   a   LATERALFLOW%ITS     aP  @   a   LATERALFLOW%ITE     ¡P  @   a   LATERALFLOW%JTS     áP  @   a   LATERALFLOW%JTE     !Q  @   a   LATERALFLOW%KTS     aQ  @   a   LATERALFLOW%KTE    ¡Q  ÷       UPDATEWTD     R  @   a   UPDATEWTD%NSOIL    ØR  Ä   a   UPDATEWTD%DZS     S  ä   a   UPDATEWTD%ZSOIL     T  Ä   a   UPDATEWTD%SMCEQ !   DU  @   a   UPDATEWTD%SMCMAX !   U  @   a   UPDATEWTD%SMCWLT !   ÄU  @   a   UPDATEWTD%PSISAT    V  @   a   UPDATEWTD%BEXP    DV  @   a   UPDATEWTD%ILOC    V  @   a   UPDATEWTD%JLOC #   ÄV  @   a   UPDATEWTD%TOTWATER    W  @   a   UPDATEWTD%WTD    DW  Ä   a   UPDATEWTD%SMC    X  Ä   a   UPDATEWTD%SH2O !   ÌX  @   a   UPDATEWTD%SMCWTD "   Y  @   a   UPDATEWTD%QSPRING 