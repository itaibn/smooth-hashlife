#include "shlife.h"

// GMP doesn't have a maximum or minimum function for mpz_t, and they'll be
// needed later
void my_mpz_max(mpz_t rop, const mpz_t op0, const mpz_t op1) {
    if (mpz_cmp(op0, op1) > 0) {
        mpz_set(rop, op0);
    } else {
        mpz_set(rop, op1);
    }
}

void my_mpz_min(mpz_t rop, const mpz_t op0, const mpz_t op1) {
    if (mpz_cmp(op0, op1) > 0) {
        mpz_set(rop, op1);
    } else {
        mpz_set(rop, op0);
    }
}


//// TOOLS FOR CALCULATING THE HASH FUNCTION
// The intended hash function is as follows: Given a block b of size 2^nx2^n, it
// sums modulo hashprime for all 0 =< i,j < 2^n the value xmul^i*ymul^j if the
// coordinate (i,j) is dead in b and 2*xmul^i*ymul^j if (i, j) is alive in b.
// This way simply multiplying the hash by xmul^n*ymul^m "translates" the block
// by (n, m), which allows calculating the hashes of subblocks more easily.

// Currently very little thought was put in choosing hashprime, xmul, & ymul.
// hashprime was chosen as a large prime that fits in 32 bits, and xmul & ymul
// were chosen so as not to have any obvious linear relationships which might
// cause collisions. It would be good to have that both xmul and ymul are
// primitive roots but this hasn't been checked.

// Note: hashprime must be less than 2^31 for code to work. See comment in
// hash_node().

int
init_hash_cache() {
    int i, j, x, y;
    int p, tmp;
    unsigned long yhash, xyhash, hash;

    mpz_init_set_ui(hashprime_mpz, hashprime);

    yhash = 1;
    for (y=0; y<LEAFSIZE; y++) {
        xyhash = yhash;
        for (x=0; x<LEAFSIZE; x++) {
            point_hash_value[y*LEAFSIZE + x] = xyhash;
            xyhash = xyhash * xmul % hashprime;
        }
        yhash = yhash * ymul % hashprime;
    }

    for (p=0; p < 1 << (LEAFSIZE*LEAFSIZE); p++) {
        hash = 0;
        for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
            if (p & (1 << i)) {
                hash = (hash + 2*point_hash_value[i]) % hashprime;
            } else {
                hash = (hash + 1*point_hash_value[i]) % hashprime;
            }
        }
        leaf_hash_cache[p] = hash;
    }

    int size;
    unsigned long *table, point_value;
    for (x=1; x<=LEAFSIZE; x++) {
    for (y=1; y<=LEAFSIZE; y++) {
        size = 1 << (x*y);
        table = (unsigned long *) malloc(size * sizeof(unsigned long));
        for (p=0; p < size; p++) {
            hash = 0;
            for (i=0; i<x; i++) {
            for (j=0; j<y; j++) {
                point_value = point_hash_value[i + LEAFSIZE*j];
                if (p & (1 << (i+x*j))) {
                    hash = (hash + 2*point_value) % hashprime;
                } else {
                    hash = (hash + 1*point_value) % hashprime;
                }
                tmp = tmp >> 1;
            }}
            assert(hash < hashprime);
            table[p] = hash;
        }
        rectangle_hash_cache_size[(x-1) + LEAFSIZE*(y-1)] = size;
        rectangle_hash_cache[(x-1) + LEAFSIZE*(y-1)] = table;
    }
    }

    hash = point_hash_value[LEAFSIZE-1] * xmul % hashprime;
    for (i=0; i < 256; i++) {
        xmul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }

    hash = point_hash_value[LEAFSIZE*(LEAFSIZE-1)] * ymul % hashprime;
    for (i=0; i<256; i++) {
        ymul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }

    return 0;
}

// Note: d is the depth of the subnodes, not the combined node.
unsigned long
hash_node(unsigned long hnw, unsigned long hne, unsigned long hsw, unsigned long
        hse, depth_t d) {
    if (d >= 256) {
        fprintf(stderr, "This implementation currently does not supported sizes"
            "larger than 2^257\n");
        exit(1)
        //return NULL;
    }

    uint64_t xmul_d, ymul_d, xymul_d;
    xmul_d = xmul_cache[d];
    ymul_d = ymul_cache[d];
    xymul_d = xmul_d * ymul_d % hashprime;
    // Since hashprime is less than 2^31, all the terms in the sum are less than
    // 2^62, guaranteeing that there is no overflow.
    return (hnw + xmul_d*hne + ymul_d*hsw + xymul_d*hse) % hashprime;
}

// Here comes the most complicated function among the hashing utilities.
// Compute the hash of a rectangular subblock with northwest corner (x0, y0) and
// southeast corner (x1, y1). The rectangle is truncated if either (x0, y0) or
// (x1, y1) extend past the ends of base. If the parameter 'adjust' is set
// nonzero, the hash is calculated as if (x0, y0) is the origin.
unsigned long
hash_rectangle(block *base, const mpz_t ix0, const mpz_t ix1, const mpz_t iy0,
        const mpz_t iy1, int adjust) {
    //assert (mpz_cmp(ix0, ix1) =< 0 && mpz_cmp(iy0, y1) =< 0)

    unsigned long hash;
    mpz_t tmp, blocksize, zero, x0, x1, y0, y1;

    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(blocksize, LEAFSIZE);
    mpz_mul_2exp(blocksize, blocksize, base->depth);
    mpz_inits(tmp, x0, x1, y0, y1, NULL);
    my_mpz_max(x0, ix0, zero);
    my_mpz_max(y0, iy0, zero);
    my_mpz_min(x1, ix1, blocksize);
    my_mpz_min(y1, iy1, blocksize);
    //assert (mpz_sgn(x0) >= 0 && mpz_cmp(x1, blocksize) =< 0
    //    && mpz_sgn(y0) >= 0 && mpz_cmp(y1, blocksize) =< 0);
    
    if (mpz_cmp(x0, x1) >= 0 || mpz_cmp(y0, y1) >= 0) {
        hash = 0;
        goto end;
    }
    if (mpz_sgn(x0) == 0 && mpz_cmp(x1, blocksize) == 0 && mpz_sgn(y0) == 0 &&
        mpz_cmp(y1, blocksize) == 0) {
        hash = base->hash;
        goto end;
    }

    if (base->tag == LEAF_B) {
        unsigned long x0l, x1l, y0l, y1l;
        size_t size;
        unsigned long *table;
        leaf xmask, mask, pos, row, rect;
        int i;
        x0l = mpz_get_ui(x0);
        x1l = mpz_get_ui(x1);
        y0l = mpz_get_ui(y0);
        y1l = mpz_get_ui(y1);
        i = (x1l-x0l-1) + LEAFSIZE*(y1l-y0l-1);
        assert (0 <= i && i < LEAFSIZE*LEAFSIZE);
        table = rectangle_hash_cache[i];
        size = rectangle_hash_cache_size[i];
        pos = base->content.b_l;
        rect = 0;
        xmask = ((1 << (x1l - x0l)) - 1) << x0l;
        for (i = y0l; i < y1l; i++) {
            mask = xmask << (i*LEAFSIZE);
            row = (pos & mask) >> (i*LEAFSIZE+x0l);
            rect |= row << ((i-y0l)*(x1l-x0l));
        }
        assert(rect < 1<<((x1l-x0l)*(y1l-y0l)));
        assert(rect < size);
        assert(size == 1<<((x1l-x0l)*(y1l-y0l)));
        hash = point_hash_value[x0l + LEAFSIZE*y0l] * table[rect] % hashprime;
        goto end;
    } else if (base->tag == CONTAIN_B) {
        subblock s = base->content.b_c;
        block *super;
        mpz_t superx0, superx1, supery0, supery1, x_adj, y_adj;
        uint64_t xy_adj;

        super = s.superblock;
        mpz_inits(superx0, superx1, supery0, supery1, x_adj, y_adj, NULL);
        mpz_add(superx0, x0, s.x);
        mpz_add(superx1, x1, s.x);
        mpz_add(supery0, y0, s.y);
        mpz_add(supery1, y1, s.y);

        hash = hash_rectangle(super, superx0, superx1, supery0, supery1, 0);

        mpz_set_ui(x_adj, xmul);
        mpz_neg(tmp, x0);
        mpz_powm(x_adj, x_adj, tmp, hashprime_mpz);
        mpz_init_set_ui(y_adj, ymul);
        mpz_neg(tmp, y0);
        mpz_powm(y_adj, y_adj, tmp, hashprime_mpz);
        xy_adj = mpz_get_ui(x_adj) * mpz_get_ui(y_adj) % hashprime;
        hash = (unsigned long) ((xy_adj * (uint64_t) hash) % hashprime);

        mpz_clears(superx0, superx1, supery0, supery1, NULL);
        goto end;
    } else if (base->tag != NODE_B) {
        fprintf(stderr, "CONTAIN_B not supported as input to hash_rectangle()");
    }

    node n = base->content.b_n;
    unsigned long hnw, hne, hsw, hse; 
    mpz_t halfblock;

    mpz_t shiftx0, shiftx1, shifty0, shifty1;
    mpz_inits(halfblock, shiftx0, shiftx1, shifty0, shifty1, NULL);
    mpz_tdiv_q_2exp(halfblock, blocksize, 1);
    mpz_sub(shiftx0, x0, halfblock);
    mpz_sub(shiftx1, x1, halfblock);
    mpz_sub(shifty0, y0, halfblock);
    mpz_sub(shifty1, y1, halfblock);
    hnw = hash_rectangle(NW(n), x0, x1, y0, y1, 0);
    hne = hash_rectangle(NE(n), shiftx0, shiftx1, y0, y1, 0);
    hsw = hash_rectangle(SW(n), x0, x1, shifty0, shifty1, 0);
    hse = hash_rectangle(SE(n), shiftx0, shiftx1, shifty0, shifty1, 0);
    
    hash = hash_node(hnw, hne, hsw, hse, base->depth-1);
    mpz_clears(halfblock, shiftx0, shiftx1, shifty0, shifty1, NULL);

end:
    if (adjust) {
        mpz_t x_adj, y_adj;
        uint64_t xy_adj;
        mpz_init_set_ui(x_adj, xmul);
        mpz_neg(tmp, x0);
        mpz_powm(x_adj, x_adj, tmp, hashprime_mpz);
        mpz_init_set_ui(y_adj, ymul);
        mpz_neg(tmp, y0);
        mpz_powm(y_adj, y_adj, tmp, hashprime_mpz);
        xy_adj = mpz_get_ui(x_adj) * mpz_get_ui(y_adj) % hashprime;
        hash = (unsigned long) ((xy_adj * (uint64_t) hash) % hashprime);
        mpz_clears(x_adj, y_adj, NULL);
    }
    mpz_clears(tmp, blocksize, zero, x0, x1, y0, y1, NULL);
    return hash;
}

//// BASIC BLOCK CREATION FUNCTIONS

// Allocate a block with a given hash from the hash table
block * new_block(unsigned long hash);

block *
mkblock_leaf(leaf l) {
    unsigned long hash = leaf_hash_cache[l];
    block *b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = LEAF_B;
    b->content.b_l = l;
    b->depth = 0;
    return b;
}

// Combine four blocks into a node block
block *
mkblock_node(block *nw, block *ne, block *sw, block *se) {
    unsigned long hash;
    block *b;
    depth_t d;

    if (nw == NULL || ne == NULL || sw == NULL || se == NULL) {
        TRACE("null case\n");
        return NULL;
    }
    d = nw->depth;
    if (ne->depth != d || sw->depth != d || se->depth != d) {
        TRACE("bad sizes\n");
        return NULL;
    }

    node n = {{{nw, ne}, {sw, se}}};
    hash = hash_node(nw->hash, ne->hash, sw->hash, se->hash, d);
    b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = NODE_B;
    b->content.b_n = n;
    b->depth = d+1;
    return b;
}

block *block_index(block *b, int y, int x);

// TODO: Seperate this into two different functions, one for diff=1, one for
// diff>1 (possibly on diff=2 necessary). Thus far this procedure is a mixup
// between block creation and more complicated block processing.
block *
mkblock_contain(block *superblock, mpz_t x, mpz_t y, depth_t diff) {
    //TRACE("mkc ss %d d %d x %Zd y %Zd\n", (int) LGLENGTH(superblock), diff, x,
    //    y);
    //gmp_printf("mkc s(%lu %p) x %Zd y %Zd dd %lu\n", (unsigned long)
    //    superblock->depth, superblock, x, y, (unsigned long) diff);
    //display(superblock, stdout);
    //printf("\n");
    assert(superblock);
    assert(mpz_sgn(x) >= 0);
    assert(mpz_sgn(y) >= 0);
    assert(diff >= 0);
    assert(diff <= superblock->depth);

    depth_t d;
    d = superblock->depth;

    //depth_t dd = LGLENGTH(superblock);
    //assert(dd >= diff);
    if (diff == 0) {
        assert(mpz_sgn(x) == 0 && mpz_sgn(y) == 0);
        return superblock;
    }

    // Verify (x,y) not too far south or east.
    mpz_t rjust, supersize;
    mpz_inits(rjust, supersize, NULL);

    mpz_init_set_ui(supersize, 1);
    mpz_mul_2exp(supersize, supersize, LGLENGTH(superblock));

    mpz_set_ui(rjust, 1);
    mpz_mul_2exp(rjust, rjust, LGLENGTH(superblock) - diff);
    mpz_add(rjust, rjust, x);
    assert (mpz_cmp(rjust, supersize) <= 0);
    //TRACE("xrjust %Zd ", rjust);

    mpz_set_ui(rjust, 1);
    mpz_mul_2exp(rjust, rjust, LGLENGTH(superblock) - diff);
    mpz_add(rjust, rjust, y);
    assert (mpz_cmp(rjust, supersize) <= 0);
    //TRACE("yrjust %Zd\n", rjust);

    mpz_clears(rjust, supersize, NULL);

    if (superblock->tag == CONTAIN_B) {
        mpz_t nx, ny;
        block *res;
        mpz_init(nx); mpz_init(ny);
        mpz_add(nx, x, superblock->content.b_c.x);
        mpz_add(ny, y, superblock->content.b_c.y);
        res = mkblock_contain(superblock->content.b_c.superblock, nx, ny,
            diff+1);
        mpz_clear(nx); mpz_clear(ny);
        assert(res->depth > 0 || res->tag == LEAF_B);
        return res;
    }
    
    if (diff > 1) {
        mpz_t tmp, nx, ny;
        block *subsuperblock, *res;
        unsigned long x_approx, y_approx;
        depth_t logsize = LGLENGTH(superblock);

        mpz_init(tmp);
        mpz_init(nx); mpz_init(ny);
        mpz_tdiv_q_2exp(tmp, x, logsize-2);
        x_approx = mpz_get_ui(tmp);
        mpz_mul_2exp(tmp, tmp, logsize-2);
        mpz_sub(nx, x, tmp);
        mpz_tdiv_q_2exp(tmp, y, logsize-2);
        y_approx = mpz_get_ui(tmp);
        mpz_mul_2exp(tmp, tmp, logsize-2);
        mpz_sub(ny, y, tmp);
        mpz_clear(tmp);

        assert(x_approx < 3 && y_approx < 3);
        subsuperblock = block_index(superblock, (int) y_approx, (int) x_approx);
        
        res = mkblock_contain(subsuperblock, nx, ny, diff-1);
        mpz_clear(nx); mpz_clear(ny);
        assert(res->depth > 0 || res->tag == LEAF_B);
        return res;
    }

    assert(diff == 1);
    assert(superblock->tag == NODE_B);
    
    if (d == 1) {
        //TRACE("mkleaf case of mkcontain\n");
        //TRACE("super ");
        //if (DEBUG) display_raw(superblock, stdout);
        //TRACE(" x %Zd y %Zd d %d\n", x, y, (int) diff);
        unsigned long xi, yi;
        xi = mpz_get_ui(x);
        yi = mpz_get_ui(y);
        assert(xi <= 2 && yi <= 2);
        return block_index(superblock, (int) xi, (int) yi);
    }

    unsigned long hash;
    mpz_t size, x_east, y_south;
    block *b;

    assert(d >= 0);
    mpz_init_set_ui(size, LEAFSIZE);
    mpz_mul_2exp(size, size, d-diff);
    mpz_init(x_east);
    mpz_add(x_east, x, size);
    mpz_init(y_south);
    mpz_add(y_south, y, size);

    mpz_mul_2exp(size, size, diff);
    //gmp_printf("- d %lu size %Zd x_east %Zd y_south %Zd\n", (unsigned long) d,
    //    size, x_east, y_south);
    assert(mpz_sgn(x) >= 0 && mpz_cmp(x_east, size) <= 0 &&
        mpz_sgn(y) >= 0 && mpz_cmp(y_south, size) <= 0
        && superblock);

    hash = hash_rectangle(superblock, x, x_east, y, y_south, 1);
    mpz_clears(size, x_east, y_south, NULL);

    b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = CONTAIN_B;
    b->depth = d - diff;
    b->content.b_c.superblock = superblock;
    mpz_init_set(b->content.b_c.x, x);
    mpz_init_set(b->content.b_c.y, y);
    // ... b->content.b_c.t ...
    assert(b->depth > 0 || b->tag == LEAF_B);
    return b;
}

// Given a 2^nx2^n block b, return the 2^(n-1)x2^(n-1) subblock with northwest
// corner at (jx*2^(n-2), iy*2^(n-2)). For example:
//
//      N
//  +-+---+-+
//  | |   | |
//  | | r | |
//  | |   | |
// W| +---+ |E    r = block_index(A, 0, 1)
//  |   A   |
//  |       |
//  |       |
//  +-------+
//      S
block *
block_index(block *b, int iy, int jx) {
    //printf("bindex b(%lu %lu) %d %d\n", (unsigned long) b->depth, b->hash, i,
    //    j);
    //display(b, stdout);
    //printf("\n");

    assert(0 <= iy && iy <= 2 && 0 <= jx && jx <= 2);

    if (b == NULL) {return NULL;}
    if (b->tag == LEAF_B) {return NULL;}

    if (b->tag == CONTAIN_B) {
        //printf("CONTAIN_B\n");
        mpz_t x, y, halfblock, shift;
        block *res;

        mpz_inits(x, y, halfblock, shift, NULL);
        mpz_set_ui(halfblock, 1);
        mpz_mul_2exp(halfblock, halfblock, LGLENGTH(b)-2);
        mpz_mul_ui(shift, halfblock, jx);
        mpz_add(x, b->content.b_c.x, shift);
        mpz_mul_ui(shift, halfblock, iy);
        mpz_add(y, b->content.b_c.y, shift);
        res = mkblock_contain(b->content.b_c.superblock, x, y, 2);
        mpz_clears(x, y, halfblock, shift, NULL);
        
        return res;
    } else if (b->tag != NODE_B) {
        fprintf(stderr, "Invalid tag %d for block at %p (hash %lu)\n", b->tag,
            b, b->hash);
    }

    assert(b->tag == NODE_B);
    node b_no = b->content.b_n;
    
    if (((iy | jx) & 1) == 0) {
        //printf("Easy corner\n");
        // This isn't just for optimization; block_index(b, i, j) where either i
        // or j is odd recursively calls the case where both i and j are even.
        return CORNER(b_no, iy/2, jx/2);
    } else if (b->depth == 1) {
    // The next two cases have similar structures. The block b is determined to
    // be a node. The subblock is computed by computing its (p, q)'th corner as
    // (p, q) varies from (0,0) to (1,1). i0 is the higher order bit of p+i and
    // i1 is the lower order bit, and similarly (j0,j1) are the bits of q+j.
    //
    // The two cases are for when the b has depth one and the result is a leaf
    // and when b has depth greater than one and the result is a node. Perhaps
    // the proper way to do this is to have a 'subleaf' block which consists of
    // a quarter of a leaf as well as a general method for combining blocks into
    // nodes as well as subleafs into leafs.
        //printf("depth 1 node\n");

        leaf res, bit;
        int p, i0, i1, q, j0, j1;
        res = 0;
        for (p=0; p<2; p++) {
        for (q=0; q<2; q++) {
            i0 = ((p + iy) & 2) >> 1;
            i1 = (p + iy) & 1;
            j0 = ((q + jx) & 2) >> 1;
            j1 = (q + jx) & 1;
            assert(CORNER(b_no, i0, j0)->tag == LEAF_B);
            bit = 1 & (CORNER(b_no, i0, j0)->content.b_l >> (LEAFSIZE*i1 + j1));
            res |= bit << (LEAFSIZE*p + q);
        }}
        //printf("l %x\n", res);
        return mkblock_leaf(res);
    } else {
        //printf("depth >1 node\n");
        node n;
        block *tmpb;
        int p, i0, i1, q, j0, j1;
        for (p=0; p<2; p++) {
        for (q=0; q<2; q++) {
            i0 = ((p + iy) & 2) >> 1;
            i1 = (p + iy) & 1;
            j0 = ((q + jx) & 2) >> 1;
            j1 = (q + jx) & 1;
            tmpb = CORNER(b_no, i0, j0);
            CORNER(n, p, q) = block_index(tmpb, 2*i1, 2*j1);
        }}
        return mkblock_node(NW(n), NE(n), SW(n), SE(n));
    }
}

//// CA COMPUTATION PROPER

int
copy_inner_pattern(inner_pattern *a, inner_pattern *b) {
    a->depth_diff = b->depth_diff;
    a->pattern = b->pattern;
    mpz_init_set(a->y, b->y);
    mpz_init_set(a->x, b->x);
}

// This code was written from scratch and not tested so surely fails badly. It
// also could use some better organizing, I think.
int
add_foci(block *b) {
    assert(b);
    if (b->nfocus >= 0) {
        return b->nfocus;
    }

    block *tmp;
    struct inner_pattern foci[58008];
    int i, j, k, count;
    count = 0;

    if (b->depth < 3) {
        if (b->depth < 2) {
            b->foci = NULL
            return (b->nfocus = 0);
        }

        mpz_t iz, jz;
        mpz_inits(iz, jz, NULL);
        for (j=2; j<7; j++) {
        for (i=2; i<7; i++) {
            mpz_set_ui(jz, j-1);
            mpz_set_ui(iz, i-1);
            foci[count].pattern = mkblock_contain(b, iz, jz, 2);
            foci[count].depth_diff = 2;
            mpz_inits(foci[count].y, foci[count].x, NULL);
            mpz_add_ui(foci[count].y, jz, 1);
            mpz_add_ui(foci[count].x, iz, 1);
            count++;
        }}
        mpz_clears(iz, jz, NULL);
    } else if (b->tag == CONTAIN_B) {
        mpz_t mnorth, msouth, mwest, meast, lhs, hsize;
        mpz_inits(mnorth, msouth, meast, mwest, lhs, hsize, NULL);

        mpz_set_ui(mnorth, 1);
        mpz_mul_2exp(hsize, mnorth, LGLENGTH(b));
        mpz_mul_2exp(mnorth, mnorth, LGLENGTH(b)-2); // Am I sure it should be 2
        mpz_set(mwest, mnorth);                      // here?
        mpz_sub(msouth, hsize, mnorth);
        mpz_sub(meast, hsize, mnorth);               // PS. This is a horrible
        mpz_add(mnorth, mnorth, b->content.b_c.y);   // way to format comments.
        mpz_add(msouth, msouth, b->content.b_c.y);
        mpz_add(mwest, mwest, b->content.b_c.x);
        mpz_add(meast, meast, b->content.b_c.x);
        mpz_tdiv_q_2exp(hsize, hsize, 1);

        for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            tmp = block_index(b->content.b_c.superblock, i, j);
            if (add_foci(tmp) < 0) {return -1; /*MPZ NOT DEALLOCATED!!!*/}
            for (k=0; k<tmp->nfocus; k++) {
                int cond;
                /*
                fx, fy = tmp->foci[k].{x,y} + {i,j}
                x, y = b->content.b_c.{x,y}
                cond = fx + j*length(b) >= x + length(b)/4
                    && fx + j*length(b) <= x + length(b)*3/4
                    && fy + i*length(b) >= y + length(b)/4
                    && fy + i*length(b) <= y + length(b)*3/4
                */
                mpz_mul_ui(lhs, hsize, i);
                mpz_add(lhs, lhs, tmp->foci[k].y);
                cond = (mpz_cmp (lhs, mnorth) >= 0) && (mpz_cmp (lhs, msouth)
                    <= 0);
                mpz_mul_ui(lhs, hsize, j);
                mpz_add(lhs, lhs, tmp->foci[k].x);
                cond = cond && (mpz_cmp (lhs, mwest) >= 0) && (mpz_cmp (lhs,
                    meast) <= 0);
                if (cond) {
//                    foci[count++] = tmp->foci[k];
                    foci[count].pattern = tmp->foci[k].pattern;
                    foci[count].depth_diff = tmp->foci[k].depth_diff;
                    mpz_init(foci[count].y);
                    mpz_mul_ui(lhs, hsize, j);
                    mpz_add(foci[count].y, tmp->foci[k].y, lhs);
                    mpz_sub(foci[count].y, foci[count].y, b->content.b_c.y);
                    mpz_init(foci[count].x);
                    mpz_mul_ui(lhs, hsize, i);
                    mpz_add(foci[count].x, tmp->foci[k].x, lhs);
                    mpz_sub(foci[count].x, foci[count].x, b->content.b_c.x);
                    count++;
                }
            }
        }}
        mpz_clears(mnorth, msouth, meast, mwest, lhs, hsize, NULL);

        b->nfocus = count;
        b->foci = (struct inner_pattern *) malloc(count * sizeof(struct
            inner_pattern));
        for (i=0; i<count; i++) {
            b->foci[i] = foci[i];
        }
    } else if (b->tag == NODE_B) {
        for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            tmp = block_index(b, i, j);
            if (tmp->nfocus < 0) {
                if (add_foci(tmp) < 0) {
                    return -1;
                }
                assert(tmp->nfocus >= 0);
            }
            for (k=0; k<tmp->nfocus; k++) {
                copy_inner_pattern(&foci[count+k], &tmp->foci[k]);
            }
            count += k;
        }}

        for (i=0; i<count; i++) {
            mpz_sub(TTY, foci[i].y, TTSIZE);
            mpz_sub(TTX, foci[i].x, TTSIZE);
            if ((mpz_sgn(TTY) >= 0) && (mpz_sgn(TTX) >= 0)) {
                foci[i].pattern = mkblock_contain(b, TTX, TTY, 2);
            } else {
                foci[i] = {-1, NULL};
            }
        }

        int cond;
        for (i=0; i<count; i++) {
            cond = (mpz_cmp(foci[i].y, TTLMARGIN) >= 0) && (mpz_cmp(foci[i].y,
                TTRMARGIN) <= 0);
            cond = cond && (mpz_cmp(foci[i].x, TTLMARGIN) >= 0) &&
                (mpz_cmp(foci[i].x, TTRMARGIN) <= 0);
            if (cond) {
                for (k=0; k<count; k++) {
                    if (k==i || foci[k]->pattern == NULL) break;
                    mpz_sub(TTDIFFY, foci[i].y, foci[k].y);
                    mpz_sub(TTDiFFX, foci[i].x, foci[k].x);
                    mpz_abs(TTDIFFY, TTDIFFY);
                    mpz_abs(TTDIFFX, TTDIFFX);
                    if ((mpz_cmp(TTDIFFY, TTDIFFBND) <= 0) && (mpz_cmp(TTDIFFX,
                            TTDIFFBND) <= 0)) {
                        cond = cond && (foci[k].pattern->hash <
                            foci[i].pattern->hash)
                    }
                }
            }
            if (!cond) {
                
            }
        }
    }
}

void
init_result() {
    int result33[65536];
    int pos, i, bitcnt, tmp;
    for (i=0; i<512; i++) {
        pos = (i&0b111) | ((i&0b111000) << 1) | ((i&0b111000000) << 2);
        tmp = i & ~0x20;
        while (tmp > 0) {
            bitcnt += tmp & 1;
            tmp = tmp >> 1;
        }
        bitcnt = __builtin_popcount(pos & ~0x20);
        int res;
        if (i & 0x10) {
            res = bitcnt == 2 || bitcnt == 3;
        } else {
            res = bitcnt == 3;
        }
        if (res) {
            result33[pos] = 1;
        } else {
            result33[pos] = 0;
        }
    }
    for (i=0; i<65536; i++) {
        result[i] = result33[i&0x0777] |
            (result33[(i&0x0eee) >> 1] << 1) |
            (result33[(i&0x7770) >> 4] << 2) |
            (result33[(i&0xeee0) >> 5] << 3);
    }
}

block *evolve(block *x);

block *
half_evolve(block *b, int i, int j) {
    //TRACE("he b ");
    //if (DEBUG) display_raw(b, stdout);
    //TRACE(" i %d j %d\n", i, j);

    assert(b);
    assert(0 <= i && i < 2 && 0 <= j && j < 2);

    node n;
    block *tmp;
    int x,y;

    for (x=0; x<2; x++) {
    for (y=0; y<2; y++) {
        tmp = evolve(block_index(b, i+x, j+y));
        assert(tmp != NULL);
        CORNER(n, x, y) = tmp;
    }}

    return mkblock_node(NW(n), NE(n), SW(n), SE(n));
}

// Sorry, the code here is a bit messy and repetitive.
block *
evolve(block *x) {
    static int nevolve = 0;
    static int cache = 0;
    static int success;
    if (x->depth > 3) {
        TRACE("evolve %d %d\n", (int) x->depth, nevolve);
    }
    nevolve++;
    block *r;
    if (x == NULL) {
        fprintf(stderr, "NULL input to evolve()\n");
    }
    if (x->res) {
        TRACE("cache  %d %d\n", (int) x->depth, nevolve);
        cache++;
        return x->res;
    }
    if (x->tag == CONTAIN_B) {
        mpz_t tmp;
        block *new_sblock;
        unsigned long x_approx, y_approx;

        mpz_init(tmp);
        mpz_tdiv_q_2exp(tmp, x->content.b_c.x, LGLENGTH(x));
        x_approx = mpz_get_ui(tmp);
        mpz_tdiv_q_2exp(tmp, x->content.b_c.y, LGLENGTH(x));
        y_approx = mpz_get_ui(tmp);
        mpz_clear(tmp);

        assert((x_approx | y_approx) < 2);

        assert(x->content.b_c.superblock->tag == NODE_B);
        new_sblock = half_evolve(x->content.b_c.superblock, (int) x_approx,
            (int) y_approx);
        r = mkblock_contain(new_sblock, x->content.b_c.x, x->content.b_c.y, 1);
        goto end;
    }
    if (x->tag != NODE_B) {
        fprintf(stderr, "Only nodes have evolve() implemented\n");
        return NULL;
    }
    if (x->depth == 1) {
        // Subblocks are all leafs
        node n = x->content.b_n;
        assert(NW(n)->tag == LEAF_B);
        int unpack_x = (NW(n)->content.b_l & 3) |
                       ((NE(n)->content.b_l & 3) << 2) |
                       ((NW(n)->content.b_l & 12) << 2) |
                       ((NE(n)->content.b_l & 12) << 4) |
                       ((SW(n)->content.b_l & 3) << 8) |
                       ((SE(n)->content.b_l & 3) << 10) |
                       ((SW(n)->content.b_l & 12) << 10) |
                       ((SE(n)->content.b_l & 12) << 12);
        r = mkblock_leaf(result[unpack_x]);
    } else {
        int i, j;
        node n;
        block *test;
        assert(x->tag == NODE_B);
        for (i=0; i<2; i++) {
        for (j=0; j<2; j++) {
            CORNER(n, i, j) = evolve(half_evolve(x, i, j));
        }}
        r = mkblock_node(NW(n), NE(n), SW(n), SE(n));
        test = mkblock_node(NW(n), NE(n), SW(n), SE(n));
        if (r != test) {
            TRACE("cache error %d %d %p | %d %d %p [%p]\n", r->index, r->hash,
                r, test->index, test->hash, test, hashtable[test->index]);
        } else {
            TRACE("cache success %d\n", success);
        }
/*
        // Half-sized subblocks of x on the north, south, east, west, and
        // center:
        block *n, *s, *w, *e, *c;
        node no = x->content.b_n;
        assert(NW(no)->tag == NODE_B);
        node tmp; // Not an actual node; just a convient way to store four
                  // blocks associated with corners.
        // Recall mkblock_node(nw, ne, sw, se)
        // This part is tedious and error-prone
        n = mkblock_node(NE(NW(no)->content.b_n),
                         NW(NE(no)->content.b_n),
                         SE(NW(no)->content.b_n),
                         SW(NE(no)->content.b_n));
        s = mkblock_node(NE(SW(no)->content.b_n),
                         NW(SE(no)->content.b_n),
                         SE(SW(no)->content.b_n),
                         SW(SE(no)->content.b_n));
        w = mkblock_node(SW(NW(no)->content.b_n),
                         SE(NW(no)->content.b_n),
                         NW(SW(no)->content.b_n),
                         NE(SW(no)->content.b_n));
        e = mkblock_node(SW(NE(no)->content.b_n),
                         SE(NE(no)->content.b_n),
                         NW(SE(no)->content.b_n),
                         NE(SE(no)->content.b_n));
        c = mkblock_node(SE(NW(no)->content.b_n),
                         SW(NE(no)->content.b_n),
                         NE(SW(no)->content.b_n),
                         NW(SE(no)->content.b_n));
        n = evolve(n);
        s = evolve(s);
        w = evolve(w);
        e = evolve(e);
        c = evolve(c);
        NW(tmp) = evolve(mkblock_node(evolve(NW(no)), n, w, c));
        NE(tmp) = evolve(mkblock_node(n, evolve(NE(no)), c, e));
        SW(tmp) = evolve(mkblock_node(w, c, evolve(SW(no)), s));
        SE(tmp) = evolve(mkblock_node(c, e, s, evolve(SE(no))));
        r = mkblock_node(NW(tmp), NE(tmp), SW(tmp), SE(tmp));
*/
    }
    end:
    if (x->depth > 4) {
        TRACE("evolve %d %d %d %d\n", x->depth, nevolve, cache, success);
    }
    return (x->res = r);
}

//// READING AND WRITING BLOCKS

block *write_bit(block *b, unsigned long y, unsigned long x, char bit);

block *
read_life_105(FILE *f) {
    int c;
    char bit;

    while ((c = fgetc(f)) == '#') {
        if (c == EOF) return NULL;
        while ((c = fgetc(f)) != '\n') {
            if (c == EOF) return NULL;
        }
    }
    
    unsigned long y, x;
    block *b, *empty, *tmp;
    b = empty = mkblock_leaf(0);
    x = y = 0;
    while(1) {
        switch (c) {
            case '*':
                bit = 1;
                goto star_or_dot;
            case '.':
                bit = 0;
            star_or_dot:
                while (!(tmp = write_bit(b, y, x, bit))) {
                    b = mkblock_node(b, empty, empty, empty);
                    empty = mkblock_node(empty, empty, empty, empty);
                }
                b = tmp;
                x++;
                break;
            case '\r':
                if ((c = fgetc(f)) != '\n') {
                    return NULL;
                }
            case '\n':
                x = 0;
                y++;
                break;
            case EOF:
                return b;
                break;
            default:
                return NULL;
                break;
        }
        c = fgetc(f);
    }
}

int display_raw(block *b, FILE *f);

block *
blank_block(depth_t lglength) {
    if (lglength < LGLEAFSIZE) {
        fprintf(stderr, "too small argument to blank_block()\n");
    }

    block *blank;
    blank = mkblock_leaf(0);
    while (LGLENGTH(blank) < lglength) {
        blank = mkblock_node(blank, blank, blank, blank);
    }
    return blank;
}

// To aid a lazy implementation, a maximum length for lines
#define MAXMCLINES 1000000
block *
read_mc(FILE *f) {
    //return NULL; // Screw this. I should probably only start working on reading
    //             // .mc files after I have the hash-table code to support it.

    char c;

    // Read first line: "[M2] ..."
    while ((c = getc(f)) != '\n') {
        if (c == EOF) {
            return NULL;
        }
    } 
    // Read second line: "#..."
    while ((c = getc(f)) != '\n') {
        if (c == EOF) {
            return NULL;
        }
    } 

    block *blank_eight, **blocktable;
    blank_eight = mkblock_leaf(0);
    blank_eight = mkblock_node(blank_eight, blank_eight, blank_eight,
        blank_eight);
    blank_eight = mkblock_node(blank_eight, blank_eight, blank_eight,
        blank_eight);
    if (LEAFSIZE != 2)
        {fprintf(stderr, "read_mc() needs LEAFSIZE=2\n"); return NULL;}
    blocktable = malloc(MAXMCLINES * sizeof(block *));
    blocktable[0] = NULL;
    blocktable[1] = NULL;
    //if (DEBUG) display_raw(blank_eight);

    //int inw, ine, isw, ise, depth;
    int depth;
    int i,j;
    int index;
    node n;
    index = 1;
    //TRACE("i %d\n", index);
    while(index < MAXMCLINES) {
        /*
        if (DEBUG) {
            //TRACE("previous patterns (index=%d):\n", index);
            for (i=0; i<index; i++) {
                //TRACE("(%d)", i);
                display_raw(blocktable[i], stderr);
            }
            //TRACE("\n");
        }
        */
        block *current;
        c = getc(f);
        if (c == '$' || c == '*' || c == '.') {
            int x=0, y=0, bit;
            //TRACE("8x8\n");
            current = blank_eight;
            while (c != '\n') {
                switch(c) {
                    case '.':
                        bit = 0;
                        goto dot_or_star;
                    case '*':
                        bit = 1;
                    dot_or_star:
                        current = write_bit(current, y, x, bit);
                        x++;
                        break;
                    case '$':
                        x = 0;
                        y++;
                        break;
                    default:
                        fprintf(stderr, "Invalid formatting in Macrocell "
                            "file\n");
                        return NULL; // WARNING: Doesn't deallocate resources.
                }
                c = getc(f);
            }
        } else if ('0' <= c && c <= '9') {
            //TRACE("big block\n");
            int x[5];
            for (i=0; i<5; i++) {
                x[i] = 0;
                do {
                    x[i] = x[i]*10 + (c - '0');
                    c = getc(f);
                } while ('0' <= c && c <= '9');
                if (i<4) {
                    //TRACE("verify c=%c (%d)\n", c, c);
                    assert(c == ' ');
                    c = getc(f);
                    assert('0' <= c && c <= '9');
                }
                //TRACE("numeral %d\n", x[i]);
            }
            if (c == '\r') {c = getc(f);}
            assert(c == '\n');
            assert(x[1] < index && x[2] < index && x[3] < index && x[4] <
                index);

            depth = x[0];
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                if (x[2*i+j+1] == 0) {
                    CORNER(n, i, j) = blank_block(depth-1);
                } else if (x[2*i+j+1] < index) {
                    CORNER(n, i, j) = blocktable[x[2*i+j+1]];
                } else {
                    fprintf(stderr, "Bad formatting\n");
                    return NULL; // Ditto
                }
            }}
            /*
            NW(n) = blocktable[x[1]];
            NE(n) = blocktable[x[2]];
            SW(n) = blocktable[x[3]];
            SE(n) = blocktable[x[4]];
            */

            //TRACE("making block ");
            /*
            if (DEBUG) {
                for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    display(CORNER(n, i, j), stderr);
                }}
            }
            */
            //TRACE("\n");
            current = mkblock_node(NW(n), NE(n), SW(n), SE(n));
            assert(current && LGLENGTH(current) == depth);
        } else if (c == EOF) {
            break;
        } else {
            fprintf(stderr, "Nonsensical line in .mc");
        }

        blocktable[index] = current;
        index++;
    }

    TRACE("done\n");
    return blocktable[index-1];
}

// Useful for debugging.
int
display_raw(block *b, FILE *f) {
    if (!b) {
        fprintf(f, "NULL");
        return 1;
    }

    node n;
    block *btmp;
    switch (b->tag) {
        case LEAF_B:
            fprintf(f, "%x", b->content.b_l);
            break;
        case CONTAIN_B:
            gmp_fprintf(f, "c[%Zd %Zd]", b->content.b_c.x, b->content.b_c.y);
            display_raw(b->content.b_c.superblock, f);
            break;
        case NODE_B:
            fprintf(f, "n(");
            int i, j;
            n = b->content.b_n;
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                btmp = CORNER(n, i, j);
                display_raw(btmp, f);
            }}
            fprintf(f, ")");
            break;
        default:
            fprintf(stderr, "Invalid tag for block input for display()\n");
            return 1;
    }
    return 0;
}

// Used in read_life_105. Returns NULL when (x,y) is out of range.
block *
write_bit(block *b, unsigned long y, unsigned long x, char bit) {
    unsigned long size = 2;
    /*
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    */
    size = 2 << b->depth;

    if (x >= size) {
        //return -1;
        return NULL;
    }
    if (y >= size) {
        //return -2;
        return NULL;
    }

    if (b->tag == LEAF_B) {
        leaf mask = 1 << (2*y+x);
        leaf new_leaf = (b->content.b_l & ~mask) | (bit ? mask : 0);
        return mkblock_leaf(new_leaf);
    } else if (b->tag == NODE_B) {
        node n = b->content.b_n;
        if (y < size/2) {
            if (x < size/2) {
                NW(n) = write_bit(NW(n), y, x, bit);
            } else {
                NE(n) = write_bit(NE(n), y, x-size/2, bit);
            }
        } else {
            if (x < size/2) {
                SW(n) = write_bit(SW(n), y-size/2, x, bit);
            } else {
                SE(n) = write_bit(SE(n), y-size/2, x-size/2, bit);
            }
        }
        if (NW(n) == NULL || NE(n) == NULL || SW(n) == NULL || SE(n) == NULL) {
            return NULL;
        } else {
            return mkblock_node(NW(n), NE(n), SW(n), SE(n));
        }
    } else {
        fprintf(stderr, "No code for writing in a block neither leaf nor node");
        return NULL;
    }
}

// Display the line from (x0, y) to (x1, y) in Life 1.05 format onto file f.
// Ends with a newline if newline!=0.
int
display_line(block *b, mpz_t y, mpz_t x0, mpz_t x1, int newline, FILE *f) {
    assert(mpz_cmp(x0, x1) <= 0);
    mpz_t size, halfsize, tmpx0, tmpx1, tmpy;

    mpz_init(size); mpz_init(halfsize); mpz_init(tmpx0); mpz_init(tmpx1);
    mpz_init(tmpy);
    mpz_set_ui(size, LEAFSIZE);
    mpz_mul_2exp(size, size, b->depth);
    mpz_tdiv_q_2exp(halfsize, size, 1);

    if (mpz_sgn(y) < 0 || mpz_cmp(y, size) >= 0 || mpz_sgn(x1) <= 0 ||
        mpz_cmp(x0, size) >= 0) {
        goto end;
    }

    unsigned long yl, x0l, x1l;
    int i, j, bit;
    leaf bt;
    node nb;
    subblock sb;
    switch (b->tag) {
        case LEAF_B:
            yl = mpz_get_ui(y);
            if (mpz_sgn(x0) < 0) {
                x0l = 0;
            } else {
                x0l = mpz_get_ui(x0);
            }
            if (mpz_cmp(x1, size) > 0) {
                x1l = LEAFSIZE;
            } else {
                x1l = mpz_get_ui(x1);
            }
            bt = b->content.b_l;

            assert(x0l < x1l);
            for (i=x0l; i<x1l; i++) {
                bit = 1 & (bt >> (yl*LEAFSIZE + i));
                fputc(bit ? '*' : '.', f);
            }
            break;

        case NODE_B:
            nb = b->content.b_n;
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                mpz_mul_si(tmpy, halfsize, -i);
                mpz_add(tmpy, tmpy, y);
                mpz_mul_si(tmpx0, halfsize, -j);
                mpz_add(tmpx0, tmpx0, x0);
                mpz_mul_si(tmpx1, halfsize, -j);
                mpz_add(tmpx1, tmpx1, x1);

                display_line(CORNER(nb, i, j), tmpy, tmpx0, tmpx1, 0, f);
            }}
            break;

        case CONTAIN_B:
            sb = b->content.b_c;

            mpz_set_ui(tmpy, 0);
            my_mpz_max(tmpx0, x0, tmpy);
            my_mpz_min(tmpx1, x1, size);

            mpz_add(tmpy, y, sb.y);
            mpz_add(tmpx0, tmpx0, sb.x);
            mpz_add(tmpx1, tmpx1, sb.x);

            display_line(sb.superblock, tmpy, tmpx0, tmpx1, 0, f);
            break;

        default:
            fprintf(stderr, "Invalid block type\n");
            return 1;
    }

    end:
    if (newline) {
        fputc('\n', f);
    }

    mpz_clear(size); mpz_clear(halfsize); mpz_clear(tmpx0); mpz_clear(tmpx1);
    mpz_clear(tmpy);
    return 0;
}

// Displays a block in Life 1.05 format.
int
display(block *b, FILE *f) {
    int ret=0;
    int size;
    size = LEAFSIZE << b->depth;
    if (size == 0) {
        fprintf(stderr, "Can't display block: Too big\n");
        return 1;
    }

    //TRACE("d size %d", size);
    //display_raw(b, stdout);
    //TRACE("\n");

    fputs("#Life 1.05\n", f);

    int i;
    mpz_t y, x0, x1;
    mpz_init_set_ui(y, 0);
    mpz_init_set_ui(x0, 0);
    mpz_init_set_ui(x1, size);
    for (i=0; i<size; i++) {
        //TRACE("d l %Zd ", y);
        if (display_line(b, y, x0, x1, 1, f)) {
            fprintf(stderr, "Error display line %d of block\n", i);
            ret = 1;
            goto end;
        }
        mpz_add_ui(y, y, 1);
    }
end:
    mpz_clears(x0, x1, y, NULL);
    return ret;
}

/*
int
print_line(block *b, long y, FILE *f) {
    unsigned long size = 2;
    / *
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    * /
    size = 2 << b->depth;

    if (b->tag == LEAF_B) {
        int out = 0xf & b->content.b_l >> (2*y);
        int j;
        for (j = 0; j<2; j++) {
            fputc((out&1) ? '*' : '.', f);
            out >>= 1;
        }
    } else if (b->tag == NODE_B) {
        if (y < size/2) {
            print_line(NW(b->content.b_n), y, f);
            print_line(NE(b->content.b_n), y, f);
        } else {
            print_line(SW(b->content.b_n), y - size/2, f);
            print_line(SE(b->content.b_n), y - size/2, f);
        }
    }
}

int
display(block *b, FILE *f) {
    if (b->tag != LEAF_B && b->tag != NODE_B) {
        fprintf(stderr, "CONTAIN_B not supported\n");
        return;
    }

    unsigned long size = 2;
    / *
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    * /
    size = 2 << b->depth;
    if (size == 0) {
        fprintf(stderr, "Overflow error: code doesn't support blocks of this"
            "size\n");
        return 1;
    }

    long i;
    for (i = 0; i < size; i++) {
        print_line(b, i, f);
        fputc('\n', f);
    }
}
*/

//// HASH TABLE-RELATED FUNCTIONS
// Note: Currently this sucks, there is no hash table resizing nor garbage
// collection.

/*
unsigned int ht_size = 1000000;
block **hashtable;
*/

void
init_hashtable() {
    hashtable = (block **) malloc(ht_size * sizeof(block *));
    memset(hashtable, 0, ht_size * sizeof(block *));
}

block *
new_block(unsigned long hash) {
    unsigned long ii = 0;
    int i;
    block *b;
    do {
        ii = ii + hash;
        i = ii % ht_size;
        b = hashtable[i];
        if (b == NULL) {
            b = (block *) malloc(sizeof(block));
            b->hash = hash;
            b->res = NULL;
            b->tag = EMPTY;
            b->nfocus = -1;
            b->foci = NULL;
#if DEBUG
            b->index = i;
#endif
//            b->refcount = 0;
            // Insert command here of the form:
            //  num_entries++;
            hashtable[i] = b;
        }
    } while (b->hash != hash);
//    b->refcount++;
    return b;
}

/*
// Decrease the refrence count of a block.
int
unref(block *b) {
    if (--(b->refcount)) {
        return 1;
    } else {
        
*/

//// MAIN

int
main(int argc, char **argv) {
    //TRACE("EMPTY %d LEAF_B %d NODE_B %d CONTAIN_B %d\n", EMPTY, LEAF_B, NODE_B,
    //    CONTAIN_B);

    init_hashtable();
    init_hash_cache();
    init_result();

    if (argc != 7) {
        fprintf(stderr, "shlife pattern outfile x0 x1 y0 y1\n");
        exit(1);
    }

    FILE *f;
    block *b;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(1);
    }
    //b = read_life_105(f);
    b = read_mc(f);
    if (b == NULL) {
        fprintf(stderr, "Badly formatted input\n");
        exit(1);
    } else if (b->depth == 0) {
        fprintf(stderr, "Input pattern must be larger than leaf\n");
        // Not exiting intentionally.
    }
    fclose(f);

    mpz_t x[4];
    int i, j, err=0;
    for (i=0; i<4; i++) {
        err |= mpz_init_set_str(x[i], argv[i+3], 10);
    }
    if (err) {
        fprintf(stderr, "Arguments x0, x1, y0, y1 must be numbers\n");
        exit(1);
    }

    f = fopen(argv[2], "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening output file\n");
        exit(1);
    }
    printf("evolving\n");
    b = evolve(b);
    printf("done\n");
    if (b == NULL) {
        fprintf(stderr, "Error evolving life pattern\n");
        exit(1);
    }
    display(b, f);
    fclose(f);

    exit(0);

/*
    display(b, stdout);
    
    node n;
    mpz_t tmpx, tmpy, shift;
    mpz_inits(tmpx, tmpy, shift, NULL);
    mpz_set_ui(shift, LEAFSIZE);
    mpz_mul_2exp(shift, shift, b->depth-2);
    for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
        mpz_set(tmpx, x[0]);
        mpz_addmul_ui(tmpx, shift, i);
        mpz_set(tmpy, x[1]);
        mpz_addmul_ui(tmpy, shift, j);
        CORNER(n, j, i) = mkblock_contain(b, tmpx, tmpy, 2);
    }}
    mpz_clears(tmpx, tmpy, shift, NULL);
    b = mkblock_node(NW(n), NE(n), SW(n), SE(n));
    //b = mkblock_contain(b, x[0], x[1], 1);

    //display(evolve(b), stdout);
    display(b, stdout);
    //if (b->tag == NODE_B) display(b->content.b_n.nw, stdout);
    printf("%lu\n", b->hash);
    printf("%lu\n", hash_rectangle(b, x[0], x[1], x[2], x[3], 1));
    block *c;
    for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
        printf("Subblock %d %d:\n", i, j);
        c = block_index(b, i, j);
        display(c, stdout);
    }}
    for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
        printf("Half-evolve %d %d:\n", i, j);
        c = half_evolve(b, i, j);
        display(c, stdout);
    }}
    printf("Result:\n");
    c = evolve(b);
    display(c, stdout);
    exit(0);
*/
}
