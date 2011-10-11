/**
 * \file message.h
 * Public parts of General Messaging Utility.
 *//*
 *  Created : 26 Feb 2010
 *  Author  : Hazel Marsden
 *
 *  Copyright (C) 2010 European Bioinformatics Institute
 *
 *  This file is part of the AYB base calling software.
 *
 *  AYB is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AYB is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AYB.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MESSAGE_H_
#define MESSAGE_H_

#include <stdbool.h>
#include "utility.h"            // for CSTRING

/** Quick debug output of string and float. For temporary use. */
#define DEBUGF(S,F) message(E_GENERIC_SF, MSG_DEBUG, S, F)


/**
 * Message type determines output text. Suffix indicates type of arguments required.
 * Enumeration is ordered according to the parameters required to enable bulk testing.
 * Dummy enumerations (E_END_) inserted between parameter types also used in testing.
 * See also MSG_TEXT defined in message.c.
 */
typedef enum MsgTypeT {E_DEFAULTBLOCK,
                       E_NOBLOCKS,
                       E_NOPATTERN,
                       E_BAD_ITER,
                       E_BAD_MATRIXIN,
                       E_BAD_RUNOPT,
                       E_END_NONE,
                       E_NOMEM_S,
                       E_MSG_LEVEL_S,
                       E_INPUT_DIR_S,
                       E_OUTPUT_DIR_S,
                       E_INPUT_FOUND_S,
                       E_NOPATTERN_FILE_S,
                       E_BAD_INPUT_S,
                       E_DATABLOCK_FAIL_S,
                       E_MATRIX_FAIL_S,
                       E_NOCREATE_S,
                       E_ZERO_LAMBDA_S,
                       E_END_S,
                       E_BAD_DIR_SS,
                       E_NOCREATE_DIR_SS,
                       E_CREATED_DIR_SS,
                       E_BAD_INPUT_SS,
                       E_NOINPUT_SS,
                       E_OPEN_FAIL_SS,
                       E_LANETILE_SS,
                       E_OPT_SELECT_SS,
                       E_BAD_TXT_SS,
                       E_BAD_NUM_SS,
                       E_END_SS,
                       E_BAD_CHAR_SC,
                       E_END_SC,
                       E_PATTERN_MATCH_SD,
                       E_OPT_SELECT_SD,
                       E_END_SD,
                       E_MATRIXINIT_SDD,
                       E_END_SDD,
                       E_OPT_SELECT_SE,
                       E_END_SE,
                       E_BAD_NUC_C,
                       E_END_C,
                       E_PROCESS_FAIL_D,
                       E_CYCLESIZE_D,
                       E_END_D,
                       E_CYCLESIZE_DD,
                       E_TILESIZE_DD,
                       E_INIT_FAIL_DD,
                       E_PROCESS_DD,
                       E_END_DD,
                       E_GENERIC_SS,
                       E_GENERIC_SD,
                       E_GENERIC_SU,
                       E_GENERIC_SX,
                       E_GENERIC_SF,
                       E_GENERIC_SFFFF,
                       E_DEBUG_SSD_S,
                       E_DEBUG_SSD_SS,
                       E_DEBUG_SSD_SD} MSGTYPE;

/**
 * Message severity determines level at which message appears.
 * See also MSG_SEV_TEXT defined in message.c.
 */
typedef enum MsgSeverityT {MSG_NONE, MSG_FATAL, MSG_ERR, MSG_INFO, MSG_WARN, MSG_DEBUG, MSG_NUM} MSGSEV;


/* function prototypes */

MSGSEV get_message_level(void);
CSTRING get_message_path(void);
int message(MSGTYPE type, MSGSEV sev, ...);

bool set_message_level(const char *levelstr);
void set_message_path(const CSTRING path);

bool startup_message(void);
void tidyup_message(void);

#endif /* MESSAGE_H_ */
