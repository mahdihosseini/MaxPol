function [i]=normal(i)
i=(i-min(i(:)))./range(i(:));