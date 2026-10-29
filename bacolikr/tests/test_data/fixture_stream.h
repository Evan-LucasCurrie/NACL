#ifndef BAC_FIXTURE_STREAM_H
#define BAC_FIXTURE_STREAM_H
#include <stddef.h>
#include <stdio.h>
#include <string.h>
/* Reconstruct embedded reference records in a standard C temporary stream. */
static FILE *bac_fixture_stream(const char *const *chunks, size_t count) {
    FILE *stream = tmpfile();
    if (!stream) return NULL;
    for (size_t i = 0; i < count; i++) {
        if (fputs(chunks[i], stream) == EOF) {
            fclose(stream);
            return NULL;
        }
    }
    rewind(stream);
    return stream;
}
#endif
