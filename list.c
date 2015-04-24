#include "shlife.h"

struct list_struct;

typedef struct list_struct {
    struct inner_pattern elem;
    struct list_struct *next;
} blocklist;

blocklist *insert(struct inner_pattern x, blocklist *y) {
    blocklist *pair = malloc(sizeof(blocklist));
    pair->elem = x;
    pair->next = y;
    return pair;
}

void delete(blocklist *l) {
    if (l) {
        delete(l->next);
        free(l);
    }
}

int length(blocklist *l) {
    if (!l) return 0;
    int n = length(l->next);
    return n+1;
}
