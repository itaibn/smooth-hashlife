-- Alternate implementation of the hash function in haskell
module Main where

hashprime, xmul, ymul :: Integer
hashprime = 1000000007
xmul = 2331
ymul = 121212121

xm = flip mod hashprime . (*) xmul
ym = flip mod hashprime . (*) ymul

hash :: [[Integer]] -> Integer
hash = foldr yjoin 0 . map (foldr xjoin 0)
    where
        yjoin a b = (a + b*ymul) `mod` hashprime
        xjoin a b = (a + b*xmul) `mod` hashprime

main = interact ((++"\n").show.hash.map (map toInt).lines)
    where
        toInt '.' = 1
        toInt '*' = 2
