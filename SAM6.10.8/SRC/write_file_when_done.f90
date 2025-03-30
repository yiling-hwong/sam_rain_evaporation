     
subroutine write_file_when_done
	
implicit none
    open(54,file='./nstop_has_been_reached',status='new',form='formatted')
    write(54,*) ' salut '
    close(unit=54)

end
