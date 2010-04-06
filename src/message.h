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


/** Message type determines output text. See also MSG_TEXT defined in message.c. */
typedef enum MsgTypeT {E_NOCYCLES, E_NOPATTERN, E_PATTERN_SELECT_S, E_NOMEM_S, E_NOINPUT_S, E_INPUTFOUND_S,
                       E_NODIR_SS, E_CYCLE_SELECT_D, E_CYCLESIZE_DD,
                       E_GENERIC_SS, E_GENERIC_SD, E_GENERIC_SU, E_GENERIC_SF, E_GENERIC_SFFFF,
                       E_DEBUG_SSD_S, E_DEBUG_SSD_SS, E_DEBUG_SSD_SD} MSGTYPE;

/** Message severity determines level at which message appears. See also MSG_SEV_TEXT defined in message.c. */
typedef enum MsgSeverityT {MSG_NONE, MSG_FATAL, MSG_ERR, MSG_WARN, MSG_INFO, MSG_DEBUG, MSG_NUM} MSGSEV;


/* function prototypes */

int message(MSGTYPE type, MSGSEV sev, ...);
bool set_message_level(const char *levelstr);
void set_message_path(const char *path);

void startup_message(const char *prefix);
void tidyup_message();
void tidy_message();

#endif /* MESSAGE_H_ */
