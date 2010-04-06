/*
 *  File    : message.h
 *  Created : 26 Feb 2010
 *  Author  : Hazel Marsden
 *  Purpose : Header containing public parts of General Messaging Utility
 *            Includes print macros
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


/* message type determines output text; see also MSG_TEXT defined in message.c */
enum MsgTypeT {E_INT_D, E_FLOAT_F, E_STRING_S, E_INT_DD};
extern const char *MSG_TEXT[];

/* message severity; see also MSG_SEV_TEXT defined in message.c */
enum MsgSeverityT {MSG_NONE, MSG_FATAL, MSG_ERR, MSG_WARN, MSG_INFO, MSG_NUM};
extern const char *MSG_SEV_TEXT[];
/* selected level of messages */
extern int Msg_Level;


/*
 * macros to define message functions as a print to stderr
 *
 * this is being done this way to allow parameters of varying types but avoid having
 * fprintf scattered around the code. Could be re-implemented if I work out how;
 * messages with 1 or 2 parameters, extend if required
 */
#define MSG1(typ, sev, p1)  {if (Msg_Level >= sev) \
                               { fprintf(stderr, MSG_TEXT[typ], MSG_SEV_TEXT[sev], p1); } }
#define MSG2(typ, sev, p1, p2) {if (Msg_Level >= sev) \
                                 { fprintf(stderr, MSG_TEXT[typ], MSG_SEV_TEXT[sev], p1, p2); } }


/* function prototypes */
bool set_message_level(const char *levelstr);
void set_message_path(const char *path);
void start_message(const char *prefix);
void tidyup_message();
void tidy_message();

#endif /* MESSAGE_H_ */
