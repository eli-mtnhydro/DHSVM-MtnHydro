
#ifndef DHSVM_ERROR_H
#define DHSVM_ERROR_H

extern char errorstr[];
void ReportError(char *ErrorString, int ErrorCode);
void ReportWarning(char *ErrorString, int ErrorCode);

#endif
