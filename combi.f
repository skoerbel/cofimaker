        module combi
        implicit none

        contains

c---------------------------------------------------------------------

        integer function facult(n)
        implicit none
        integer, intent(in) :: n
        integer i

        facult=1
        do i=2,n
          facult=facult*i
        end do

        end function facult

c---------------------------------------------------------------------

        function facult_large(n)
        implicit none
        integer, intent(in) :: n
        integer(selected_int_kind(15)) :: facult_large
        integer i

        facult_large=1
        do i=2,n
          facult_large=facult_large*i
        end do

        end function facult_large

c---------------------------------------------------------------------

        integer function dblfactorial(n)
        implicit none
        integer, intent(in) :: n
        integer i

        dblfactorial=1
        do i=3,n,2
          dblfactorial=dblfactorial*i
        end do

        end function

c---------------------------------------------------------------------

        integer function n_over_k(n,k)
        use defs, only : error_stop
        implicit none
        integer, intent(in) :: n,k
        !
        if (n.lt.0.or.k.lt.0)                                             &
     &           call error_stop("Non-negative integer required")
        !if (n.lt.k) call error_stop("n must be >= k")
        if (n.lt.k) n_over_k=0
        if (n.eq.k) n_over_k=1
        if (n.gt.k) then
          n_over_k=facult(n)/(facult(k)*facult(n-k))
        end if
        !
        end function n_over_k

c---------------------------------------------------------------------

        integer function n_over_k_large(n,k)
        use defs, only : error_stop
        implicit none
        integer, intent(in) :: n,k
        double precision facult_n,facult_k,facult_n_k
        integer(selected_int_kind(15)) n_large,k_large
        !
        if (n.lt.0.or.k.lt.0)                                             &
     &           call error_stop("Non-negative integer required")
        !if (n.lt.k) call error_stop("n must be >= k")
        if (n.lt.k) n_over_k_large=0
        if (n==k) n_over_k_large=1
        if (n.gt.k) then
          !n_over_k=facult(n)/(facult(k)*facult(n-k))
          facult_n=facult_large(n)
          facult_k=facult_large(k)
          facult_n_k=facult_large(n-k)
          n_over_k_large=int(facult_n/(facult_k*facult_n_k))
        end if
        !
        end function n_over_k_large

c---------------------------------------------------------------------

        end module combi
