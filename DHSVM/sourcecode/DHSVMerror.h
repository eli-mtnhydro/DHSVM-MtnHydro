#ifndef DHSVM_ERROR_H
#define DHSVM_ERROR_H

extern char errorstr[];
void ReportError(char *ErrorString, int ErrorCode);
void ReportErrorCode(int Code, const char *Format, ...);
void ReportWarning(char *ErrorString, int ErrorCode);
void Warning(const char *Format, ...);
void SetVerbose(int Verbose);
void Message(const char *Format, ...);

#endif
