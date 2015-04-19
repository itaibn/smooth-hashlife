-- Alternate implementation of the hash function in haskell

-- A haskell model for certain aspects of shlife to validate the code with a
-- less optimized and simpler implementation
module Model where

import Prelude hiding (length, (!!))
import Data.List (genericLength, genericIndex)

-- Integers for everything
type I = Integer
length = genericLength
(!!) = genericIndex

type Block = [[I]]

readBoard :: String -> Block
readBoard = map (map toInt) . lines
    where
        toInt '.' = 1
        toInt '*' = 2
        toInt _ = error "Badly formed life position"

size :: Block -> I
size = length

leafsize :: I
leafsize = 2

hashprime, xmul, ymul :: I
hashprime = 1000000007
xmul = 2331
ymul = 121212121

xm = flip mod hashprime . (*) xmul
ym = flip mod hashprime . (*) ymul

hash :: Block -> I
hash = foldr yjoin 0 . map (foldr xjoin 0)
    where
        yjoin a b = (a + b*ymul) `mod` hashprime
        xjoin a b = (a + b*xmul) `mod` hashprime

depthToSize :: I -> I
depthToSize n = leafsize * 2^n

subblock :: I -> Block -> (I,I) -> Block
subblock m b (i,j) = let
    s = 2^m
    n = i-2^(m-1)
    w = j-2^(m-1)
  in
    --[b !! x !! y | x <- [n..n+s], y <- [w..w+s]]
    [[b !! y !! x | x <- [w..w+s-1]] | y <- [n..n+s-1]]

-- (Apologies for the vague pronouns and specifiers in the followng.)
--
-- An i-focal subblock of a block must have size 2^i*leafsize and have the
-- 2^(i+1)*leafsize subblock which contains it at the center is entirely
-- enclosed by the superblock. Any subblock satisfying these conditions for i=0
-- is 0-focal. For i>0, it must satisfy that the 2^(i-1)*leafsize subblock at
-- the center of this subblock is (i-1)-focal, and that there is no (i-1)-focal
-- block whose enclosing 2^i*leafsize block has a larger hash than this in the
-- enclosing 2^(i+1)*leafsize subblock of this.
--
-- The center of a block is the north-west-most point of the south west
-- quadrant. focal i b returns the centers of the i-focal subblocks of the block
-- b.
focal :: I -> Block -> [(I, I)]
{-
focal n b = if not (checkSize (n+2) b) then error "Invalid board"
    else if n==0 then
        [b!!i!!j | i <- [], j <- []]
-}
focal 0 b = let s = size b in
    [(i,j) | i <- [1..s-1], j <- [1..s-1]] -- ad hoc simplification of the interval
focal i b = let
    siz = size b
    base = focal (i-1) b
    posCond (n,m) = n >= depthToSize i && n < siz-depthToSize i && m >=
        depthToSize i && m < siz-depthToSize i
    posFiltered = filter posCond base
    isNearby (n,m) (t,u) = max (abs (n-t)) (abs (m-u)) <= depthToSize i `div` 2
    challengers foc = filter (isNearby foc) $ filter (/=foc) base
    rating pt = hash (subblock i b pt)
    isFocal foc = all (<rating foc) (map rating (challengers foc))
    in filter isFocal posFiltered

{-
main = interact ((++"\n").show.hash.map (map toInt).lines)
    where
        toInt '.' = 1
        toInt '*' = 2
-}
