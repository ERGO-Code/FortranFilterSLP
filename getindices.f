      INTEGER FUNCTION RESULTINDEX(SPEC,INGR,ITEM,NUMCOMPS,NUMITEMS)
     &                                
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : RESULTINDEX                          |
*       | Description: Returns the index for Spec,NT/Rm,Item|
*       | Revision   : September 2003                       |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------------
*  Calculate the index  *
*------------------------
*
        RESULTINDEX=(SPEC-1)*NUMCOMPS*NUMITEMS
     &              +(INGR-1)*NUMITEMS
     &              +ITEM
        RETURN
        END




      INTEGER FUNCTION RELINDEX(GIG,ITEM,NUMNTS,NUMRMS,NUMPMX,
     &                            NUMSPS,NRATIOS,NUMGIG,NUMITEMS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : RELINDEX                             |
*       | Description: Returns the index for REL, Item      |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------------
*  Calculate the index  *
*------------------------
*
        OFFSET  =(NUMPMX+NUMSPS)*(NUMNTS+NUMRMS+NUMPMX)*8
     &          +(NUMRMS+NUMPMX)*8
     &          + NRATIOS*8
     &          + NUMGIG*8
*
        RELINDEX=(GIG-1)*NUMITEMS
     &          +ITEM
     &          +OFFSET
        RETURN
        END






      INTEGER FUNCTION RATIOINDEX(RATIO,ITEM,NUMNTS,NUMRMS,NUMPMX,
     &                              NUMSPS,NUMITEMS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : RATIOINDEX                           |
*       | Description: Returns the index for Ratio,Item     |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------------
*  Calculate the index  *
*------------------------
*
        OFFSET    =(NUMPMX+NUMSPS)*(NUMNTS+NUMRMS+NUMPMX)*8
     &            +(NUMRMS+NUMPMX)*8
*
        RATIOINDEX=(RATIO-1)*NUMITEMS
     &             +ITEM
     &             +OFFSET
        RETURN
        END






      INTEGER FUNCTION GIGINDEX(GIG,ITEM,NUMNTS,NUMRMS,NUMPMX,
     &                            NUMSPS,NRATIOS,NUMITEMS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : GIGINDEX                             |
*       | Description: Returns the index for GIG, Item      |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------------
*  Calculate the index  *
*------------------------
*
        OFFSET  =(NUMPMX+NUMSPS)*(NUMNTS+NUMRMS+NUMPMX)*8
     &          +(NUMRMS+NUMPMX)*8
     &          + NRATIOS*8
*
        GIGINDEX=(GIG-1)*NUMITEMS
     &          +ITEM
     &          +OFFSET
        RETURN
        END







      INTEGER FUNCTION INFINDEX(REL,ITEM,NUMNTS,NUMRMS,NUMPMX,
     &                            NUMSPS,NRATIOS,NUMGIG,NUMREL,NUMITEMS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : INFINDEX                             |
*       | Description: Returns the index for INF, Item      |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------------
*  Calculate the index  *
*------------------------
*
        OFFSET  =(NUMPMX+NUMSPS)*(NUMNTS+NUMRMS+NUMPMX)*8
     &          +(NUMRMS+NUMPMX)*8
     &          + NRATIOS*8
     &          + NUMGIG*8
     &          + NUMREL*8
*
        INFINDEX=(REL-1)*NUMITEMS
     &          +ITEM
     &          +OFFSET
        RETURN
        END
